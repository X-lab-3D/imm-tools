# SPDX-License-Identifier: Apache-2.0
"""Shared input handling and PDB line utilities (pdb-tools conventions).

Every tool reads a PDB or mmCIF file (see _doc) given on the command line or on stdin, takes at most
one option written as -<value>, and writes to stdout; messages go to stderr.
"""
import io
import os
import sys

COORD = ("ATOM  ", "HETATM", "ANISOU", "SIGATM", "SIGUIJ")


# residue references in header records: (chain column, residue number columns, icode column)
REFS = {
    "HELIX ": [(19, (21, 25), 25), (31, (33, 37), 37)],
    "SHEET ": [(21, (22, 26), 26), (32, (33, 37), 37), (49, (50, 54), 54), (64, (65, 69), 69)],
    "SSBOND": [(15, (17, 21), 21), (29, (31, 35), 35)],
    "CISPEP": [(15, (17, 21), 21), (29, (31, 35), 35)],
    "LINK  ": [(21, (22, 26), 26), (51, (52, 56), 56)],
}


def _refs(line):
    """[(chain, number, icode, field)] residue references of a header record."""
    out = []
    for f in REFS.get(line[:6], []):
        c, (a, b), i = f
        if len(line.rstrip("\n")) <= a:
            continue
        try:
            num = int(line[a:b])
        except ValueError:
            continue
        out.append((line[c], num, line[i].strip() if len(line) > i else "", f))
    return out


def _set_field(line, f, chain=None, num=None, icode=None):
    c, (a, b), i = f
    line = line.rstrip("\n").ljust(i + 1)
    if chain is not None:
        line = line[:c] + chain + line[c + 1:]
    if num is not None:
        line = line[:a] + f"{num:>{b - a}d}" + (icode or " ")[:1] + line[i + 1:]
    return line.rstrip() + "\n"


def check_input(args, doc, allow_empty_option=True):
    """(option, file handle) from sys.argv[1:] following pdb-tools rules."""
    option, fh = "", sys.stdin
    if not args:
        if sys.stdin.isatty():
            sys.stderr.write(doc)
            sys.exit(1)
    elif len(args) == 1:
        if args[0] in ("-h", "--help"):
            sys.stderr.write(doc)
            sys.exit(0)
        if args[0].startswith("-"):
            option = args[0][1:]
            if sys.stdin.isatty():
                sys.stderr.write("ERROR!! No data to process!\n")
                sys.stderr.write(doc)
                sys.exit(1)
        else:
            fh = _open(args[0], doc)
    elif len(args) == 2:
        if not args[0].startswith("-"):
            sys.stderr.write(f"ERROR!! First argument is not an option: '{args[0]}'\n")
            sys.stderr.write(doc)
            sys.exit(1)
        option = args[0][1:]
        fh = _open(args[1], doc)
    else:
        sys.stderr.write(doc)
        sys.exit(1)
    if not option and not allow_empty_option:
        sys.stderr.write("ERROR!! This tool needs an option.\n")
        sys.stderr.write(doc)
        sys.exit(1)
    return option, fh


def _open(path, doc):
    if not os.path.isfile(path):
        sys.stderr.write(f"ERROR!! File not found or not readable: '{path}'\n")
        sys.stderr.write(doc)
        sys.exit(1)
    return open(path, "r")


def chain_list(option):
    """'A,B' -> {'A', 'B'}; '' -> None (all chains)."""
    ids = {c.strip() for c in option.split(",") if c.strip()}
    return ids or None


def structure(lines):
    """Bio.PDB structure from PDB lines."""
    from Bio.PDB import PDBParser
    if not any(ln.startswith(("ATOM  ", "HETATM")) for ln in lines):
        sys.stderr.write("ERROR!! No ATOM/HETATM records in the input.\n")
        sys.exit(1)
    return PDBParser(QUIET=True).get_structure("input", io.StringIO("".join(lines)))


def residue_key(res):
    """Key of a Bio.PDB residue that matches line_key() of its lines."""
    return (res.get_parent().id, res.id[1], res.id[2].strip(), res.id[0] != " ")


def line_key(line, is_het):
    return (line[21], int(line[22:26]), line[26].strip(), is_het)


def apply_numbering(lines, new_ids):
    """Rewrite residue numbers in coordinate records (and TER) of `lines`.

    new_ids: {residue_key: (number, insertion_code)} for every renumbered residue.
    Residues of a renumbered chain that are not in new_ids (waters, ligands, ...) keep their
    number unless it collides with a new one; then they move after the chain's highest new
    number.  Residue references in HELIX, SHEET, SSBOND, LINK and CISPEP records are
    renumbered too (SEQRES and other header records are left unchanged)."""
    chains = {k[0] for k in new_ids}
    used = {}
    for k, v in new_ids.items():
        used.setdefault(k[0], set()).add(v)
    # residues that must move
    moved, nxt = {}, {c: max(n for n, _ in used[c]) + 1 for c in chains}
    last_het = False
    for ln in lines:
        rec = ln[:6]
        if rec in ("ATOM  ", "HETATM"):
            last_het = rec == "HETATM"
        if rec in COORD and ln[21] in chains:
            k = line_key(ln, last_het)
            if k not in new_ids and k not in moved:
                cur = (k[1], k[2])
                if cur in used[k[0]]:
                    moved[k] = (nxt[k[0]], "")
                    nxt[k[0]] += 1
                else:
                    moved[k] = cur
                used[k[0]].add(moved[k])
    # header records (HELIX, SHEET, SSBOND, LINK, CISPEP) follow the new numbering, so
    # viewers that read them (e.g. PyMOL) show the right secondary structure
    by_res = {}
    for k, v in list(moved.items()) + list(new_ids.items()):
        by_res[(k[0], k[1], k[2])] = v
    order = {}                                   # residues per chain, file order (1st model)
    for ln in lines:
        if ln.startswith("ENDMDL"):
            break
        if ln[:6] in ("ATOM  ", "HETATM") and ln[21] in chains:
            r = (ln[21], int(ln[22:26]), ln[26].strip())
            lst = order.setdefault(ln[21], [])
            if not lst or lst[-1] != r:
                lst.append(r)
    out, last_het, last_new = [], False, None
    for ln in lines:
        rec = ln[:6]
        if rec in ("HELIX ", "SHEET "):
            out.extend(_ss_record(ln, by_res, order))
            continue
        if rec in REFS:
            for ch, num, ic, f in _refs(ln):
                v = by_res.get((ch, num, ic))
                if v is not None:
                    ln = _set_field(ln, f, num=v[0], icode=v[1])
            out.append(ln)
            continue
        if rec in ("ATOM  ", "HETATM"):
            last_het = rec == "HETATM"
        if rec in COORD and ln[21] in chains:
            k = line_key(ln, last_het)
            num, ic = new_ids.get(k) or moved[k]
            if not -999 <= num <= 9999:
                raise ValueError(f"residue number {num} does not fit the PDB format")
            ln = f"{ln[:22]}{num:>4d}{(ic or ' ')[:1]}{ln[27:]}"
            last_new = (ln[21], num, ic)
        elif rec == "TER   " and len(ln) > 26 and last_new and ln[21] == last_new[0]:
            ln = f"{ln[:22]}{last_new[1]:>4d}{(last_new[2] or ' ')[:1]}{ln[27:]}"
        out.append(ln)
    return out


def ss_runs(lst, by_res, i0, i1):
    """Split the secondary-structure range lst[i0..i1] (old residue keys (chain, num, icode)
    in chain order) into runs of new ids that a viewer reads correctly: each run ascends in
    (number, insertion code) order and its interval holds no other residue of the chain
    (PyMOL takes every residue between start and end in that order). Returns [[(num, ic)]]."""
    import bisect
    key = lambda v: (v[0], v[1] or "")
    allkeys = sorted(key(by_res.get(r, (r[1], r[2]))) for r in lst)
    new = [by_res.get(r, (r[1], r[2])) for r in lst[i0:i1 + 1]]
    runs, cur = [], [new[0]]
    for v in new[1:]:
        lo, hi = key(cur[0]), key(v)
        inside = bisect.bisect_right(allkeys, hi) - bisect.bisect_left(allkeys, lo)
        if hi > key(cur[-1]) and inside == len(cur) + 1:
            cur.append(v)
        else:
            runs.append(cur)
            cur = [v]
    runs.append(cur)
    return runs


def _ss_record(line, by_res, order):
    """HELIX / SHEET record in the new numbering.

    Viewers such as PyMOL assign a HELIX / SHEET range to the residues between start and end
    in (number, insertion code) order. IMGT numbering has descending insertion codes (C-LIKE
    1H..1A, 85G..85A, CDR3 112G..112A), so one range in the old numbering can become a
    range PyMOL cannot read (e.g. 1001F..1001B). Such ranges are split into runs that ascend."""
    refs = _refs(line)
    if len(refs) < 2 or refs[0][0] != refs[1][0] or refs[0][0] not in order:
        for ch, num, ic, f in refs:
            v = by_res.get((ch, num, ic))
            if v is not None:
                line = _set_field(line, f, num=v[0], icode=v[1])
        return [line]
    ch = refs[0][0]
    lst = order[ch]
    try:
        i0, i1 = lst.index(refs[0][:3]), lst.index(refs[1][:3])
    except ValueError:
        i0 = i1 = None
    if i0 is None or i1 < i0:
        for c, num, ic, f in refs:
            v = by_res.get((c, num, ic))
            if v is not None:
                line = _set_field(line, f, num=v[0], icode=v[1])
        return [line]
    runs = ss_runs(lst, by_res, i0, i1)
    out = []
    for k, run in enumerate(runs):
        ln = _set_field(line, refs[0][3], num=run[0][0], icode=run[0][1])
        ln = _set_field(ln, refs[1][3], num=run[-1][0], icode=run[-1][1])
        if line[:6] == "HELIX " and len(ln.rstrip("\n")) >= 76:
            ln = ln[:71] + f"{len(run):>5d}" + ln[76:]
        if line[:6] == "SHEET ":
            if k == 0:                           # registration stays with the first part
                for c, num, ic, f in refs[2:]:
                    v = by_res.get((c, num, ic))
                    if v is not None:
                        ln = _set_field(ln, f, num=v[0], icode=v[1])
            else:
                ln = ln.rstrip("\n")[:38].rstrip() + "\n"
        out.append(ln)
    return out


def write(lines):
    """Write a string or a list of lines to stdout."""
    try:
        sys.stdout.write(lines if isinstance(lines, str) else "".join(lines))
        sys.stdout.flush()
    except IOError:            # e.g. piped into head
        pass


CHAIN_RECORDS = COORD + ("TER   ",)


def select_chains(lines, keep):
    """Lines of the chains in `keep` (coordinates, TER, CONECT between kept atoms); all
    other records (header, MODEL/ENDMDL, END) are kept."""
    kept_serials, out = set(), []
    for ln in lines:
        rec = ln[:6]
        if rec in CHAIN_RECORDS:
            if len(ln) > 21 and ln[21] in keep:
                out.append(ln)
                if rec in ("ATOM  ", "HETATM"):
                    kept_serials.add(ln[6:11].strip())
            continue
        if rec == "CONECT":
            ser = [ln[i:i + 5].strip() for i in range(6, min(len(ln.rstrip("\n")), 31), 5)]
            if all(s in kept_serials for s in ser if s):
                out.append(ln)
            continue
        if rec in REFS:
            if all(r[0] in keep for r in _refs(ln)):
                out.append(ln)
            continue
        out.append(ln)
    return out


def rename_chains(lines, mapping):
    """Change chain identifiers (column 22) by mapping {old: new} in chain records."""
    out = []
    for ln in lines:
        if ln[:6] in CHAIN_RECORDS and len(ln) > 21 and ln[21] in mapping:
            ln = ln[:21] + mapping[ln[21]] + ln[22:]
        elif ln[:6] in REFS:
            for ch, num, ic, f in _refs(ln):
                if ch in mapping:
                    ln = _set_field(ln, f, chain=mapping[ch])
        out.append(ln)
    return out


def resid(res):
    """'123A' style residue id of a Bio.PDB residue."""
    return f"{res.id[1]}{res.id[2].strip()}"


def name_of(fh, default="input"):
    return os.path.splitext(os.path.basename(fh.name))[0] if fh is not sys.stdin else default
