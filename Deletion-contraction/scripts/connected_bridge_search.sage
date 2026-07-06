load("isolating_utils.sage")
import itertools

n = 5
mismatches = []
tested = 0
skipped_disconnected = 0

def test_all_edges(name, maximal, F, q):
    global tested, skipped_disconnected
    Delta = complex_from_maximal(maximal)
    if not is_connected(Delta):
        skipped_disconnected += 1
        return
    F = frozenset(F)
    for e in itertools.combinations(sorted(F), 2):
        e = frozenset(e)
        tested += 1
        priv = is_private(Delta, F, e)
        iso = is_isolating(Delta, F, e)
        ok, diff, chi_con = check_theorem(Delta, F, e, n, q)
        if iso != ok:
            mismatches.append((name, q, maximal, sorted(F), sorted(e), priv, iso, ok))

# Family 1: triangle F={0,1,2} + edge B={3,4}, ALL bridge subsets (empty one now excluded by connectivity check)
A, B = (0,1,2), (3,4)
all_bridges = [(a,b) for a in A for b in B]
for q in [2, 3]:
    for k in range(0, len(all_bridges)+1):
        for combo in itertools.combinations(all_bridges, k):
            maximal = [A, B] + list(combo)
            test_all_edges("tri+edge", maximal, A, q)

# Family 2: tetrahedron F (|F|=4) + point, all bridge subsets
A2 = (0,1,2,3)
possible_bridges2 = [(v,4) for v in A2]
for q in [2,3]:
    for k in range(len(possible_bridges2)+1):
        for combo in itertools.combinations(possible_bridges2, k):
            maximal = [A2, (4,)] + list(combo)
            test_all_edges("tetra+point", maximal, A2, q)

# Family 3: edge + triangle F, all bridge subsets
B3 = (2,3,4)
possible_bridges3 = [(a,b) for a in (0,1) for b in B3]
for q in [2,3]:
    for k in range(len(possible_bridges3)+1):
        for combo in itertools.combinations(possible_bridges3, k):
            maximal = [(0,1), B3] + list(combo)
            test_all_edges("edge+tri(F=tri)", maximal, B3, q)

# Family 4: shared vertex (always connected already)
for q in [2,3]:
    maximal = [(0,1,2),(2,3,4)]
    test_all_edges("shared-vertex F=A", maximal, (0,1,2), q)
    test_all_edges("shared-vertex F=B", maximal, (2,3,4), q)
    maximal2 = [(0,1,2),(2,3,4),(0,3)]
    test_all_edges("shared-vertex+bridge F=A", maximal2, (0,1,2), q)
    test_all_edges("shared-vertex+bridge F=B", maximal2, (2,3,4), q)

print(f"Tested (connected only): {tested} (F,e) checks")
print(f"Skipped as disconnected: {skipped_disconnected} complex instances")
print(f"Mismatches: {len(mismatches)}")
for m in mismatches[:20]:
    print("  MISMATCH:", m)
