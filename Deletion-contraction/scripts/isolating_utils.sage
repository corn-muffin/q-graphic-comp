load("qarr.sage")

def one_skeleton_edges(Delta):
    return [tuple(sorted(G)) for G in Delta if len(G) == 2]

def connected_components(vertices, edges):
    parent = {v: v for v in vertices}
    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py
    for (a, b) in edges:
        union(a, b)
    return {v: find(v) for v in vertices}

def is_private(Delta, F, e):
    e = frozenset(e)
    for G in Delta:
        if e <= G and not (G <= F):
            return False
    return True

def is_connected(Delta):
    verts = sorted({v for G in Delta for v in G})
    if len(verts) <= 1:
        return True
    edges = one_skeleton_edges(Delta)
    comp = connected_components(verts, edges)
    return len(set(comp.values())) == 1

def is_isolating(Delta, F, e):
    """i,j (endpoints of e) are in different components of the 1-skeleton of
    Delta with F's OWN edges removed (not just Delta\F, which still retains
    F's boundary edges as separate elements of the downward closure)."""
    F = frozenset(F)
    all_edges = one_skeleton_edges(Delta)
    ext_edges = [ed for ed in all_edges if not frozenset(ed) <= F]
    verts = sorted({v for G in Delta for v in G})
    comp = connected_components(verts, ext_edges)
    i, j = sorted(e)
    return comp.get(i, i) != comp.get(j, j)

def check_theorem(Delta, F, e, n, q):
    """Returns True/False/None (None if not applicable due to degree mismatch)."""
    Ddel = deletion(Delta, F)
    Dcon, n_con = relabel_contiguous(contraction(Delta, e))
    chi_D = fast_characteristic_polynomial(Delta, n, q)
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
    diff = chi_D - chi_Ddel
    chi_con = fast_characteristic_polynomial(Dcon, n_con, q)
    predicted = -(q - 1) ** (len(F) - 1) * chi_con
    return (diff == predicted), diff, chi_con
