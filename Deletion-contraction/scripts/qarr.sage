"""
Core library for computing characteristic polynomials of the q-deformed
simplicial-complex arrangements S^q_Delta (Nian, "A generalization of
q-deformation of graphic arrangements to simplicial complexes", Def 3.1),

    S^q_Delta := union over faces J of Delta of
                 { ker(a_{i1} x_{i1} + ... + a_{ir} x_{ir}) : a_ij in F_q^* }

for J = {i1,...,ir} in Delta.

Vertex set is always range(n) = {0,1,...,n-1}.  A simplicial complex is
represented simply as a frozenset of frozensets (all of its faces, i.e.
already downward closed) -- helper `complex_from_maximal` builds this
from a list of maximal faces.
"""

from sage.all import *
from itertools import combinations, product


def complex_from_maximal(maximal_faces):
    """Downward closure of a list of maximal faces (each an iterable of ints)."""
    faces = set()
    for F in maximal_faces:
        F = frozenset(F)
        for r in range(1, len(F) + 1):
            for sub in combinations(sorted(F), r):
                faces.add(frozenset(sub))
    return frozenset(faces)


def full_simplex(n):
    """All nonempty subsets of {0,...,n-1} (the full (n-1)-simplex)."""
    return complex_from_maximal([range(n)])


def skeleton(Delta, k):
    """The (k-1)-skeleton: keep only faces of size <= k."""
    return frozenset(F for F in Delta if len(F) <= k)


def link(Delta, F):
    """link_Delta(F) = { G in Delta : G cap F = empty, G cup F in Delta }."""
    F = frozenset(F)
    return frozenset(G for G in Delta if G & F == frozenset() and (G | F) in Delta)


def deletion(Delta, F):
    """Delta \\ {F}: remove the single (maximal) face F."""
    F = frozenset(F)
    assert F in Delta
    return frozenset(G for G in Delta if G != F)


def contraction(Delta, F):
    """Delta / F: identify all vertices of F to a single new vertex '*'.
    We represent '*' by the special symbol -1 (guaranteed not a vertex
    label since vertices are 0..n-1)."""
    F = frozenset(F)
    STAR = -1
    new_faces = set()
    for G in Delta:
        if G & F:
            image = frozenset((G - F) | {STAR})
        else:
            image = frozenset(G)
        # take downward closure of the image too (safety; images of faces
        # of a simplicial complex under this quotient are automatically
        # downward closed as a family, but let's be safe)
        for r in range(1, len(image) + 1):
            for sub in combinations(sorted(image, key=str), r):
                new_faces.add(frozenset(sub))
    return frozenset(new_faces)


def relabel_contiguous(Delta):
    """Relabel vertices of Delta to 0..m-1 (needed after contraction, whose
    vertex set uses -1 for the new identified vertex).  Returns (new_Delta, m)."""
    verts = sorted({v for F in Delta for v in F}, key=lambda v: (v == -1, v))
    relabel = {v: i for i, v in enumerate(verts)}
    new_Delta = frozenset(frozenset(relabel[v] for v in F) for F in Delta)
    return new_Delta, len(verts)


def generate_normals(Delta, n, q):
    """List of normal vectors (tuples of length n over GF(q)) for all
    hyperplanes of S^q_Delta.  For a face J of size r, there are
    (q-1)^(r-1) hyperplanes: fix coefficient 1 at min(J), let the other
    r-1 coefficients range freely over GF(q)^*."""
    F = GF(q)
    nz = [F(a) for a in range(1, q)]  # F_q^*
    normals = []
    for J in Delta:
        J = sorted(J)
        i0 = J[0]
        rest = J[1:]
        if not rest:
            v = [F(0)] * n
            v[i0] = F(1)
            normals.append(tuple(v))
        else:
            for coeffs in product(nz, repeat=len(rest)):
                v = [F(0)] * n
                v[i0] = F(1)
                for idx, c in zip(rest, coeffs):
                    v[idx] = c
                normals.append(tuple(v))
    return normals


def characteristic_polynomial(Delta, n, q, verbose=False):
    """Build S^q_Delta as a Sage HyperplaneArrangement over GF(q) and return
    its characteristic polynomial as a univariate polynomial in t over ZZ."""
    normals = generate_normals(Delta, n, q)
    if verbose:
        print(f"  n={n}, q={q}: {len(normals)} hyperplanes")
    F = GF(q)
    names = tuple(f"x{i}" for i in range(n))
    H = HyperplaneArrangements(F, names)
    hyps = []
    for v in normals:
        hyps.append([F(0)] + list(v))  # Sage convention: constant term FIRST, then coefficients
    A = H(*hyps)
    chi = A.characteristic_polynomial()
    R = PolynomialRing(ZZ, 't')
    return R(chi)


def chi_of_maximal_faces(maximal_faces, n, q, verbose=False):
    Delta = complex_from_maximal(maximal_faces)
    return characteristic_polynomial(Delta, n, q, verbose=verbose)


# ---------------------------------------------------------------------------
# Fast custom characteristic-polynomial computation.
#
# Sage's built-in HyperplaneArrangement.characteristic_polynomial() builds
# the intersection poset in pure Python and does not scale past a few
# hundred hyperplanes.  Since S^q_Delta is a *central* arrangement, its
# characteristic polynomial equals that of the linear matroid represented
# by the normal vectors (Orlik-Terao).  We compute the lattice of flats
# directly with Sage's (M4RI-backed, fast) finite-field linear algebra,
# tracking for every flat the bitmask of which ground vectors it contains
# (its "content"), which lets us both (a) group ground vectors correctly
# when building the next rank of flats and (b) test containment between
# flats by a single integer AND, so the final Moebius-function sum never
# needs a fresh linear-algebra call.
# ---------------------------------------------------------------------------

def _echelon_key(rows, q, n):
    """Canonical (hashable) RREF key for the span of `rows` (list of tuples)."""
    if not rows:
        return ()
    M = matrix(GF(q), rows).echelon_form()
    return tuple(tuple(r) for r in M.rows() if not r.is_zero())


def fast_characteristic_polynomial(Delta, n, q, verbose=False):
    normals = generate_normals(Delta, n, q)
    m = len(normals)
    if verbose:
        print(f"  n={n}, q={q}: {m} hyperplanes (fast method)")
    F = GF(q)
    G = matrix(F, normals)  # m x n matrix, ground vectors as rows

    # flats[key] = (rank, content_bitmask)
    bottom_key = ()
    flats = {bottom_key: (0, 0)}
    frontier = [bottom_key]

    while frontier:
        next_frontier = []
        for key in frontier:
            rank, content = flats[key]
            if rank == n:
                continue
            B = matrix(F, list(key)) if key else matrix(F, 0, n)
            # reduce every ground vector modulo span(B)
            R = G.__copy__()
            for row in B.rows():
                piv = next(i for i in range(n) if row[i] != 0)
                col = R.column(piv)
                if not col.is_zero():
                    R = R - col.column() * row.row()
            groups = {}
            for idx in range(m):
                if (content >> idx) & 1:
                    continue  # already in this flat
                r = R.row(idx)
                if r.is_zero():
                    continue  # shouldn't happen (not in content but zero residual)
                piv = next(i for i in range(n) if r[i] != 0)
                inv = r[piv] ** (-1)
                canon = tuple((inv * r[i]) for i in range(n))
                groups.setdefault(canon, []).append(idx)
            for canon, idxs in groups.items():
                rep = normals[idxs[0]]
                new_key = _echelon_key(list(key) + [rep], q, n)
                new_content = content
                for idx in idxs:
                    new_content |= (1 << idx)
                if new_key not in flats:
                    flats[new_key] = (rank + 1, new_content)
                    next_frontier.append(new_key)
                else:
                    # same flat reached from a different parent: contents must agree
                    assert flats[new_key][1] == new_content, "content mismatch"
        frontier = next_frontier

    if verbose:
        print(f"  total flats: {len(flats)}")

    # Moebius function, bottom-up by rank, using bitmask containment.
    items = sorted(flats.items(), key=lambda kv: kv[1][0])  # sort by rank
    mu = {}
    lower_ranks = []  # list of (content, mu) for all flats of strictly smaller rank than current bucket
    current_rank = -1
    this_rank_batch = []
    R_ = PolynomialRing(ZZ, 't')
    t = R_.gen()
    chi = R_(0)
    for key, (rank, content) in items:
        if rank != current_rank:
            lower_ranks.extend(this_rank_batch)
            this_rank_batch = []
            current_rank = rank
        if content == 0:
            mu_here = 1
        else:
            s = 0
            for c_j, mu_j in lower_ranks:
                if (c_j & content) == c_j:
                    s += mu_j
            mu_here = -s
        mu[key] = mu_here
        this_rank_batch.append((content, mu_here))
        chi += mu_here * t ** (n - rank)
    return chi
