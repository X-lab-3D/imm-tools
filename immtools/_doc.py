# SPDX-License-Identifier: Apache-2.0
"""Structure documents: one interface for PDB and mmCIF input / output.

    doc = read(fh)                       # format detected from the content
    model = doc.model()                  # Bio.PDB model (first model) for the analysis
    doc = doc.renumbered(new_ids)        # {residue_key: (number, insertion code)}
    doc = doc.renamed({"A": "M"})        # chain IDs
    doc = doc.selected({"A", "B"})       # only these chains
    doc = doc.residues({residue_key})    # only these residues (coordinates)
    text = doc.text()                    # same format as the input

PDB files are edited line by line (pdb-tools style, see _io). mmCIF files are edited as
CIF categories with gemmi, so everything else in the file is kept: only the author fields
(auth_asym_id, auth_seq_id, pdbx_PDB_ins_code) change; label_* fields stay as they are.
"""
import io
import sys

from . import _io

# categories with author residue references: (category, [(chain col, number col, icode col)])
CIF_REFS = {
    "_struct_conf.": [("beg_auth_asym_id", "beg_auth_seq_id", "pdbx_beg_PDB_ins_code"),
                      ("end_auth_asym_id", "end_auth_seq_id", "pdbx_end_PDB_ins_code")],
    "_struct_sheet_range.": [("beg_auth_asym_id", "beg_auth_seq_id", "pdbx_beg_PDB_ins_code"),
                             ("end_auth_asym_id", "end_auth_seq_id", "pdbx_end_PDB_ins_code")],
    "_struct_conn.": [("ptnr1_auth_asym_id", "ptnr1_auth_seq_id", "pdbx_ptnr1_PDB_ins_code"),
                      ("ptnr2_auth_asym_id", "ptnr2_auth_seq_id", "pdbx_ptnr2_PDB_ins_code")],
    "_struct_mon_prot_cis.": [("auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code"),
                              ("pdbx_auth_asym_id_2", "pdbx_auth_seq_id_2", "pdbx_PDB_ins_code_2")],
    "_pdbx_poly_seq_scheme.": [("pdb_strand_id", "pdb_seq_num", "pdb_ins_code")],
    "_pdbx_nonpoly_scheme.": [("pdb_strand_id", "pdb_seq_num", "pdb_ins_code")],
}
SS_CATEGORIES = ("_struct_conf.", "_struct_sheet_range.")


def read(fh):
    """PdbDoc or CifDoc from an open file (format from the content)."""
    text = fh.read()
    head = text.lstrip()[:200]
    if head.startswith("data_") or "\n_atom_site." in text[:20000] or head.startswith("#\n"):
        return CifDoc(text)
    return PdbDoc(text.splitlines(keepends=True))


def _ic(v):
    """mmCIF insertion code value -> '' or the letter."""
    return "" if v in (None, False, "?", ".") else str(v)


class PdbDoc:
    fmt = "pdb"
    ext = ".pdb"

    def __init__(self, lines):
        self.lines = lines
        self._model = None

    def model(self):
        if self._model is None:
            self._model = next(_io.structure(self.lines).get_models())
        return self._model

    def renumbered(self, new_ids):
        return PdbDoc(_io.apply_numbering(self.lines, new_ids)) if new_ids else self

    def renamed(self, mapping):
        for new in mapping.values():
            if len(new) != 1:
                sys.stderr.write(f"ERROR!! Chain ID '{new}' does not fit the PDB format "
                                 "(one character); use mmCIF input for longer IDs.\n")
                sys.exit(1)
        return PdbDoc(_io.rename_chains(self.lines, mapping))

    def selected(self, chains):
        return PdbDoc(_io.select_chains(self.lines, set(chains)))

    def residues(self, keys):
        out, last_het = [], False
        for ln in self.lines:
            rec = ln[:6]
            if rec in ("ATOM  ", "HETATM"):
                last_het = rec == "HETATM"
            if rec in _io.COORD:
                if _io.line_key(ln, last_het) in keys:
                    out.append(ln)
            elif rec in ("TER   ", "CONECT") or rec in _io.REFS:
                continue
            else:
                out.append(ln)
        return PdbDoc(out)

    def text(self):
        return "".join(self.lines)


class CifDoc:
    fmt = "cif"
    ext = ".cif"

    def __init__(self, text):
        import gemmi
        self.doc = gemmi.cif.read_string(text)
        self.block = self.doc.sole_block()
        self._model = None
        if self.block.find_value("_atom_site.id") is None and \
                not self.block.get_mmcif_category("_atom_site."):
            sys.stderr.write("ERROR!! No _atom_site records in the mmCIF input.\n")
            sys.exit(1)

    # ------------------------------------------------------------------ helpers
    def _copy(self):
        return CifDoc(self.doc.as_string())

    def _cat(self, name):
        return self.block.get_mmcif_category(name)

    def _set(self, name, cat):
        self.block.set_mmcif_category(name, cat)

    def _atom_keys(self, cat):
        """Residue key of every _atom_site row."""
        return [(c, int(n), _ic(i), g == "HETATM") for c, n, i, g in
                zip(cat["auth_asym_id"], cat["auth_seq_id"], cat["pdbx_PDB_ins_code"],
                    cat["group_PDB"])]

    # ------------------------------------------------------------------ interface
    def model(self):
        if self._model is None:
            from Bio.PDB import MMCIFParser
            st = MMCIFParser(QUIET=True).get_structure("input", io.StringIO(self.doc.as_string()))
            self._model = next(st.get_models())
        return self._model

    def renumbered(self, new_ids):
        if not new_ids:
            return self
        d = self._copy()
        at = d._cat("_atom_site.")
        keys = d._atom_keys(at)
        chains = {k[0] for k in new_ids}
        used = {}
        for k, v in new_ids.items():
            used.setdefault(k[0], set()).add(v)
        moved, nxt = {}, {c: max(n for n, _ in used[c]) + 1 for c in chains}
        first_model = at["pdbx_PDB_model_num"][0] if "pdbx_PDB_model_num" in at else None
        order = {}
        for i, k in enumerate(keys):
            if k[0] in chains and k not in new_ids and k not in moved:
                cur = (k[1], k[2])
                if cur in used[k[0]]:
                    moved[k] = (nxt[k[0]], "")
                    nxt[k[0]] += 1
                else:
                    moved[k] = cur
                used[k[0]].add(moved[k])
            if k[0] in chains and (first_model is None or at["pdbx_PDB_model_num"][i] == first_model):
                lst = order.setdefault(k[0], [])
                r = k[:3]
                if not lst or lst[-1] != r:
                    lst.append(r)
        by_res = {k[:3]: v for k, v in list(moved.items()) + list(new_ids.items())}
        for i, k in enumerate(keys):
            v = new_ids.get(k) or moved.get(k)
            if v is not None:
                at["auth_seq_id"][i] = str(v[0])
                at["pdbx_PDB_ins_code"][i] = v[1] or None
        d._set("_atom_site.", at)
        # residue references in other categories; SS ranges are split where needed
        label_of = {}
        for k, la, ls in zip(keys, at["label_asym_id"], at.get("label_seq_id", [])):
            v = new_ids.get(k) or moved.get(k)
            if v is not None:
                label_of[(k[0],) + tuple(v)] = (la, ls)
        for name, fields in CIF_REFS.items():
            cat = d._cat(name)
            if not cat:
                continue
            if name in SS_CATEGORIES:
                cat = _split_ss(cat, fields, by_res, order, label_of)
            else:
                n = len(next(iter(cat.values())))
                for cc, nc, icc in fields:
                    if nc not in cat or cc not in cat:
                        continue
                    for i in range(n):
                        try:
                            key = (cat[cc][i], int(cat[nc][i]), _ic(cat.get(icc, [None] * n)[i]))
                        except (TypeError, ValueError):
                            continue
                        v = by_res.get(key)
                        if v is not None:
                            cat[nc][i] = str(v[0])
                            if icc in cat:
                                cat[icc][i] = v[1] or (None if name.startswith("_struct")
                                                       else ".")
            d._set(name, cat)
        return d

    def renamed(self, mapping):
        d = self._copy()
        at = d._cat("_atom_site.")
        at["auth_asym_id"] = [mapping.get(c, c) for c in at["auth_asym_id"]]
        d._set("_atom_site.", at)
        for name, fields in CIF_REFS.items():
            cat = d._cat(name)
            if not cat:
                continue
            for cc, _, _ in fields:
                if cc in cat:
                    cat[cc] = [mapping.get(c, c) if isinstance(c, str) else c for c in cat[cc]]
            d._set(name, cat)
        return d

    def selected(self, chains):
        d = self._copy()
        at = d._cat("_atom_site.")
        keep = [i for i, c in enumerate(at["auth_asym_id"]) if c in chains]
        d._set("_atom_site.", {k: [v[i] for i in keep] for k, v in at.items()})
        for name, fields in CIF_REFS.items():
            cat = d._cat(name)
            if not cat:
                continue
            n = len(next(iter(cat.values())))
            rows = [i for i in range(n)
                    if all(cat[cc][i] in chains for cc, _, _ in fields if cc in cat
                           and isinstance(cat[cc][i], str))]
            d._set(name, {k: [v[i] for i in rows] for k, v in cat.items()})
        return d

    def residues(self, keys):
        d = self._copy()
        at = d._cat("_atom_site.")
        ak = d._atom_keys(at)
        keep = [i for i, k in enumerate(ak) if k in keys]
        d._set("_atom_site.", {k: [v[i] for i in keep] for k, v in at.items()})
        return d

    def text(self):
        return self.doc.as_string()


def _split_ss(cat, fields, by_res, order, label_of):
    """Renumber SS ranges (_struct_conf, _struct_sheet_range); split ranges a viewer would
    misread because of descending IMGT insertion codes (see _io.ss_runs)."""
    (bc, bn, bi), (ec, en, ei) = fields
    n = len(next(iter(cat.values())))
    out = {k: [] for k in cat}
    for i in range(n):
        row = {k: v[i] for k, v in cat.items()}
        try:
            start = (row[bc], int(row[bn]), _ic(row.get(bi)))
            end = (row[ec], int(row[en]), _ic(row.get(ei)))
        except (TypeError, ValueError):
            start = end = None
        runs = None
        if start and start[0] == end[0] and start[0] in order:
            lst = order[start[0]]
            if start in lst and end in lst and lst.index(start) <= lst.index(end):
                runs = _io.ss_runs(lst, by_res, lst.index(start), lst.index(end))
        if runs is None:                    # not renumbered here: copy (mapping the ends)
            for (cc, nc, icc), key in ((fields[0], start), (fields[1], end)):
                v = by_res.get(key) if key else None
                if v is not None:
                    row[nc] = str(v[0])
                    if icc in row:
                        row[icc] = v[1] or None
            for k in out:
                out[k].append(row[k])
            continue
        for j, run in enumerate(runs):
            r = dict(row)
            r[bn], r[en] = str(run[0][0]), str(run[-1][0])
            if bi in r:
                r[bi] = run[0][1] or None
            if ei in r:
                r[ei] = run[-1][1] or None
            for side, v in (("beg", run[0]), ("end", run[-1])):
                lab = label_of.get((start[0],) + tuple(v))
                if lab and f"{side}_label_seq_id" in r:
                    r[f"{side}_label_asym_id"], r[f"{side}_label_seq_id"] = lab
            if j:
                r["id"] = f"{row['id']}_{j + 1}"
            if "pdbx_PDB_helix_length" in r:
                r["pdbx_PDB_helix_length"] = str(len(run))
            for k in out:
                out[k].append(r[k])
    return out
