# SPDX-License-Identifier: Apache-2.0
"""Group chains into complexes by chain type and contacts.

Complexes are built by ASSEMBLERS, run in order on the chains not yet used. Each assembler
returns a list of complexes (dict role -> chain id) and marks the chains it uses. To support
a new kind of complex (e.g. antibody-antigen), write an assembler and add it to ASSEMBLERS
(see CONTRIBUTING.md); the tools that use complexes (imm_complex, imm_rechain, imm_split,
imm_annotate) then handle it.
"""
from . import _chains
from ._registry import ROLE_OF_KIND, ROLES  # noqa: F401  (ROLES re-exported)

CONTACT = 8.0         # CA-CA distance (A) counted as a contact
MIN_CONTACTS = 3


def contacts(model, ids):
    """{(a, b): number of CA-CA pairs within CONTACT A} for chain pairs."""
    from Bio.PDB import NeighborSearch
    atoms = [a for c in model if c.id in ids for r in c for a in r if a.get_id() == "CA"]
    out = {}
    if not atoms:
        return out
    for a, b in NeighborSearch(atoms).search_all(CONTACT):
        ca, cb = a.get_parent().get_parent().id, b.get_parent().get_parent().id
        if ca != cb:
            k = tuple(sorted((ca, cb)))
            out[k] = out.get(k, 0) + 1
    return out


class Context:
    """Chains, their kinds and roles, contacts, and which chains are still free."""

    def __init__(self, model, ident):
        self.ident = ident
        self.role = {cid: ROLE_OF_KIND.get(c["kind"]) for cid, c in ident.items()}
        self.ids = [c.id for c in model if self.role.get(c.id)]        # chains with a role
        self.all = [c.id for c in model if ident[c.id]["kind"] != "none"]  # every protein chain
        self.C = contacts(model, set(self.all))
        self.free = set(self.all)

    def n(self, a, b):
        """Number of contacts between chains a and b."""
        return self.C.get(tuple(sorted((a, b))), 0)

    def chains(self, role):
        """Free chains with this role, in file order."""
        return [c for c in self.ids if c in self.free and self.role[c] == role]

    def free_chains(self, families):
        """Free chains of the given families ("other", "peptide", ...), in file order."""
        from ._registry import family
        return [c for c in self.all if c in self.free and family(self.ident[c]["kind"]) in families]

    def use(self, cid):
        self.free.discard(cid)

    def take(self, role, partners, cx, min_contacts=MIN_CONTACTS):
        """Add to complex cx the free chain of `role` with most contacts to `partners`."""
        best = max(self.chains(role), key=lambda c: sum(self.n(c, p) for p in partners),
                   default=None)
        if best and sum(self.n(best, p) for p in partners) >= min_contacts:
            cx[role] = best
            self.use(best)
            return best
        return None


# --------------------------------------------------------------------------- assemblers

def tcr_pmhc(ctx):
    """MHC heavy chain (or class II alpha) + beta2m (or class II beta) + peptide in the
    groove + the TCR pair touching them."""
    out = []
    for m in ctx.chains("mhc"):
        cx = {"mhc": m}
        ctx.use(m)
        ctx.take("b2m" if ctx.ident[m]["kind"] == "MHC_I" else "mhcb", [m], cx)
        groove = [cx[r] for r in ("mhc", "mhcb") if r in cx]
        ctx.take("pep", groove, cx)
        face = groove + ([cx["pep"]] if "pep" in cx else [])
        a = ctx.take("tcra", face, cx)
        if a:
            ctx.take("tcrb", [a], cx)
        else:
            b = ctx.take("tcrb", face, cx)
            if b:
                ctx.take("tcra", [b], cx)
        out.append(cx)
    return out


def lone_b2m(ctx):
    out = []
    for c in ctx.chains("b2m"):
        ctx.use(c)
        out.append({"b2m": c})
    return out


def tcr_only(ctx):
    """TCR alpha/beta pairs (and single chains) without pMHC."""
    out = []
    for a in ctx.chains("tcra"):
        ctx.use(a)
        cx = {"tcra": a}
        ctx.take("tcrb", [a], cx)
        out.append(cx)
    for b in ctx.chains("tcrb"):
        ctx.use(b)
        out.append({"tcrb": b})
    return out


# Order matters: earlier assemblers get first pick of the chains.
ASSEMBLERS = [tcr_pmhc, lone_b2m, tcr_only]


def find(model, identified=None):
    """List of complexes: dict(role -> chain id), in file order of their first chain.
    identified: {chain id: _chains.identify() result} (computed if not given)."""
    ident = identified or {c.id: _chains.identify(c) for c in model}
    ctx = Context(model, ident)
    complexes = []
    for assemble in ASSEMBLERS:
        complexes += assemble(ctx)
    order = {c: i for i, c in enumerate(ctx.all)}
    complexes = [cx for cx in complexes if cx]
    complexes.sort(key=lambda cx: min(order[c] for c in cx.values()))
    return complexes


def describe(cx):
    return " ".join(f"{r}={cx[r]}" for r in ROLES if r in cx)
