# SPDX-License-Identifier: Apache-2.0
"""
Sanity-checks the numbering of a (TCR-)pMHC structure and reports one line per check as a
tab-separated table: file, chain, type, check, status (OK / WARN / FAIL), detail.

Checks:
    numbering   the file's residue numbers vs IMGT numbering computed here (TCR V domain:
                ANARCI; MHC / beta2m: imm_mhcnum; peptide: 1..n)
    anchors     conserved residues at their IMGT positions in the file numbering:
                TCR V   Cys23, Trp41, Cys104, Phe/Trp118
                MHC I   Cys1011, Cys1074 (alpha2 disulfide), Cys2023, Trp2041, Cys2104 (alpha3)
                MHC II alpha  Cys1023, Trp1041, Cys1104;  MHC II beta  Cys11, Cys74,
                Cys1023, Trp1041, Cys1104;  beta2m  Cys23, Cys104
    confidence  MHC: lowest domain bit score, and residues numbered differently by the two
                model sets (broad, mammal-focused) - these positions are uncertain
    breaks      chain breaks (CA-CA > 4.2 A) with the estimated number of missing residues

Usage:
    imm_check [-<chain id(s)>] <pdb or mmcif file>

Example:
    imm_check 1ao7_imgt.pdb
    imm_tcrnum 1ao7.pdb | imm_mhcnum | imm_pepnum | imm_check

Exit status 0, also when checks fail (read the status column).

This program is part of imm-tools, made in the style of pdb-tools.
"""
import sys

from immtools import _io, _doc, _chains, _mhc, _tcr
from immtools._registry import ANCHORS, family


def _cmp(rs, labels):
    pairs = [(r, l) for r, l in zip(rs, labels) if l]
    bad = [(r, l) for r, l in pairs if (r.id[1], r.id[2].strip()) != (l[0], l[1])]
    n = len(pairs)
    pct = 100 * (n - len(bad)) / n if n else 0
    status = "OK" if not bad else ("WARN" if pct >= 90 else "FAIL")
    ex = ", ".join(f"{_io.resid(r)}->{l[0]}{l[1]}" for r, l in bad[:5])
    detail = f"{n - len(bad)}/{n} residues as IMGT ({pct:.1f}%)" + (f"; e.g. {ex}" if bad else "")
    return status, detail


def _anchors(rs, kind):
    from immtools._chains import sequence
    key = "V" if family(kind) in ("tcr", "antibody") else kind
    if key not in ANCHORS:
        return None
    by_num = {(r.id[1], r.id[2].strip()): r for r in rs}
    issues, ok = [], []
    for pos, aa in ANCHORS[key]:
        r = by_num.get((pos, ""))
        if r is None:
            issues.append(f"{pos} missing")
        elif sequence([r]) not in aa:
            issues.append(f"{pos} is {r.get_resname()}")
        else:
            ok.append(f"{sequence([r])}{pos}")
    status = "OK" if not issues else "WARN"
    return status, ("; ".join(issues) + (" | " if issues and ok else "") + " ".join(ok))


def _confidence(c):
    rs = c["residues"]
    labels, info = c["mhc"]
    worst = min(info["domains"], key=lambda d: d[3])
    sets = list(_mhc.DEFAULT_DIRS) if not (_mhc.MODEL_DIRS or _mhc.MODEL_DIR) else []
    diff = []
    if len(sets) > 1:
        saved = (_mhc.MODEL_DIR, _mhc.MODEL_DIRS)
        per = []
        try:
            for d in sets:
                _mhc.MODEL_DIRS, _mhc.MODEL_DIR = [d], None
                _mhc._cache.clear()
                try:
                    per.append(_mhc.number_residues(rs)[0])
                except ValueError:
                    per.append(None)
        finally:
            _mhc.MODEL_DIR, _mhc.MODEL_DIRS = saved
            _mhc._cache.clear()
        if all(per):
            diff = [_io.resid(r) for r, a, b in zip(rs, *per) if a != b]
    status = "OK" if not diff and worst[3] >= 100 else "WARN"
    detail = f"lowest domain score {worst[0]} {worst[3]:.0f} bits"
    if diff:
        detail += f"; {len(diff)} residues numbered differently by the two model sets: " + \
                  ",".join(diff[:15]) + (" ..." if len(diff) > 15 else "")
    return status, detail


def _breaks(rs):
    gaps = _mhc.chain_gaps([r for r in rs if "CA" in r])
    ca = [r for r in rs if "CA" in r]
    if not gaps:
        return "OK", "no chain breaks"
    items = [f"{_io.resid(ca[i - 1])}/{_io.resid(ca[i])} (~{n} missing)" for i, n in sorted(gaps.items())]
    return "WARN", f"{len(gaps)} break(s): " + ", ".join(items)


COLS = ("file", "chain", "type", "check", "status", "detail")


def run(fhandle, option, name):
    model = _doc.read(fhandle).model()
    yield "\t".join(COLS) + "\n"
    for r in rows(model, _io.chain_list(option), name):
        yield "\t".join(r[k] for k in COLS) + "\n"


def rows(model, chains=None, name=""):
    """One dict per check (keys: COLS)."""
    for ch in model:
        if chains and ch.id not in chains:
            continue
        c = _chains.identify(ch)
        if c["kind"] == "none":
            continue
        rs, kind = c["residues"], c["kind"]
        found = []
        if c["tcr"]:
            found.append(("numbering",) + _cmp(rs, c["tcr"][0]))
        elif c["mhc"]:
            found.append(("numbering",) + _cmp(rs, c["mhc"][0]))
        elif kind == "peptide":
            found.append(("numbering",) + _cmp(rs, [(k, "") for k in range(1, len(rs) + 1)]))
        a = _anchors(rs, kind)
        if a:
            found.append(("anchors",) + a)
        if c["mhc"]:
            found.append(("confidence",) + _confidence(c))
        found.append(("breaks",) + _breaks(rs))
        for check, status, detail in found:
            yield dict(file=name, chain=ch.id, type=kind, check=check, status=status,
                       detail=detail)


def main():
    option, fh = _io.check_input(sys.argv[1:], __doc__)
    _io.chain_list(option)
    if not _tcr.available():
        sys.stderr.write("WARNING: ANARCI not installed (pip install anarci-pyhmmer): TCR chains are "
                         "not checked.\n")
    _io.write(list(run(fh, option, _io.name_of(fh, "stdin"))))
    fh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
