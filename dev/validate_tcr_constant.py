"""TCR constant-domain numbering vs IMGT/3Dstructure-DB residue numbering (TR chains)."""
import sys, glob, gzip, io, re, collections
sys.path.insert(0, '.')
from immtools import _tcr
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.PDBData import protein_letters_3to1 as T31


def chains_of(text):
    out, cur = {}, None
    for ln in text.splitlines():
        m = re.match(r"REMARK 410 Chain ID\s+\S+_(\S)\s", ln)
        if m: cur = m.group(1)
        m = re.match(r"REMARK 410 IMGT chain description\s+(\S+)", ln)
        if m and cur: out[cur] = m.group(1); cur = None
    return out


tot = collections.Counter(); ok = collections.Counter(); ex = collections.Counter(); seen = set()
for f in sorted(glob.glob(sys.argv[1] + '/*.gz')):
    txt = gzip.open(f, 'rt').read(); desc = chains_of(txt)
    try:
        st = PDBParser(QUIET=1).get_structure('x', io.StringIO(txt))
    except Exception:
        continue
    for ch, d in desc.items():
        if not d.startswith('TR-') or ch not in st[0]:
            continue
        rs = [r for r in st[0][ch] if is_aa(r, standard=False)]
        seq = ''.join(T31.get(r.get_resname(), 'X') for r in rs)
        if seq in seen or not any(1000 < r.id[1] < 2000 for r in rs):
            continue
        seen.add(seq)
        res = _tcr.full_labels(seq, constant=True)
        if not res:
            continue
        lab, info = res
        k = d.split('-')[1]
        if any(r.id[1] > 3000 for r in rs):                # two V domains (scTv): skip
            tot["skipped scTv"] += 1
            continue
        tot[(k, 'chains')] += 1
        for r, l in zip(rs, lab):
            t = (r.id[1], r.id[2].strip())
            # IMGT/3Dstructure-DB: TRAC C domain 2001-2128 (one residue at 1001), TRBC 1001-
            if (1000 < t[0] < 2000 and k != 'ALPHA') or (2000 < t[0] < 3000 and k == 'ALPHA'):
                t = (1000 + t[0] % 1000, t[1])
                tot[k] += 1; ok[k] += t == l
                if t != l: ex[(k, f"{t[0]}{t[1]}", f"{l[0]}{l[1]}")] += 1
for k in ('ALPHA', 'BETA'):
    if tot[k]: print(f"{k}: {tot[(k, 'chains')]} unique chains, C-domain residues {ok[k]}/{tot[k]} = {100*ok[k]/tot[k]:.1f}%")
print(tot['skipped scTv'], 'chains with two V domains skipped'); print(ex.most_common(12))
