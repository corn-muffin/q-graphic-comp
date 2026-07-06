load("isolating_utils.sage")
import itertools

n = 5
mismatches = []
tested = 0

def test_all_edges(name, maximal, F, q):
    global tested
    Delta = complex_from_maximal(maximal)
    F = frozenset(F)
    for e in itertools.combinations(sorted(F), 2):
        e = frozenset(e)
        tested += 1
        priv = is_private(Delta, F, e)
        iso = is_isolating(Delta, F, e)
        ok, diff, chi_con = check_theorem(Delta, F, e, n, q)
        if iso != ok:
            mismatches.append((name, q, maximal, sorted(F), sorted(e), priv, iso, ok))

# Family 2: A = tetrahedron {0,1,2,3}, B = single vertex {4}, bridges subset of 4 possible edges to 4
# (4,) included explicitly so the arrangement stays essential even with 0 bridges)
A = (0,1,2,3)
possible_bridges = [(v,4) for v in A]
for q in [2,3]:
    for k in range(len(possible_bridges)+1):
        for combo in itertools.combinations(possible_bridges, k):
            maximal = [A, (4,)] + list(combo)
            test_all_edges("tetra+point", maximal, A, q)

# Family 3: A = edge {0,1}, B = triangle {2,3,4}, bridges subset of 6 possible; test F=B
B = (2,3,4)
possible_bridges3 = [(a,b) for a in (0,1) for b in B]
for q in [2,3]:
    for k in range(len(possible_bridges3)+1):
        for combo in itertools.combinations(possible_bridges3, k):
            maximal = [(0,1), B] + list(combo)
            test_all_edges("edge+tri(F=tri)", maximal, B, q)

# Family 4: A={0,1,2}, B={2,3,4} sharing vertex 2 (no extra bridges), test F=A and F=B
for q in [2,3]:
    maximal = [(0,1,2),(2,3,4)]
    test_all_edges("shared-vertex F=A", maximal, (0,1,2), q)
    test_all_edges("shared-vertex F=B", maximal, (2,3,4), q)
    # with an extra bridge between the "free" ends, e.g. 0-3
    maximal2 = [(0,1,2),(2,3,4),(0,3)]
    test_all_edges("shared-vertex+bridge F=A", maximal2, (0,1,2), q)
    test_all_edges("shared-vertex+bridge F=B", maximal2, (2,3,4), q)

print(f"Tested {tested} (complex,edge) checks. Mismatches: {len(mismatches)}")
for m in mismatches[:30]:
    print("  MISMATCH:", m)
