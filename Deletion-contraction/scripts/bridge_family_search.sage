load("isolating_utils.sage")
import itertools

n = 5

def run_case(name, A, B, bridges, q):
    """A, B partition {0,...,n-1} (as maximal faces), bridges = list of (a,b) edges."""
    maximal = [tuple(sorted(A)), tuple(sorted(B))] + [tuple(sorted(br)) for br in bridges]
    Delta = complex_from_maximal(maximal)
    F = frozenset(A)
    results = []
    for e in itertools.combinations(sorted(A), 2):
        e = frozenset(e)
        priv = is_private(Delta, F, e)
        iso = is_isolating(Delta, F, e)
        ok, diff, chi_con = check_theorem(Delta, F, e, n, q)
        results.append((e, priv, iso, ok))
    return results

mismatches = []
tested = 0

# Family 1: A = triangle {0,1,2}, B = edge {3,4}, vary which bridges (subset of the 6 possible A-B edges)
A, B = (0,1,2), (3,4)
all_bridges = [(a,b) for a in A for b in B]
for q in [2, 3]:
    for k in range(0, len(all_bridges)+1):
        for combo in itertools.combinations(all_bridges, k):
            tested += 1
            results = run_case(f"A=tri,B=edge,bridges={combo}", A, B, combo, q)
            for (e, priv, iso, ok) in results:
                if iso != ok:
                    mismatches.append(("triABedge", q, combo, sorted(e), priv, iso, ok))

print(f"Family 1 (triangle + edge, all {2**6} bridge subsets, q=2,3): tested {tested} (complex,edge) instances so far")
print(f"Mismatches so far: {len(mismatches)}")
for m in mismatches[:20]:
    print("  MISMATCH:", m)
