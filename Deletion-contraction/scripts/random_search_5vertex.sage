load("isolating_utils.sage")
import itertools, random

n = 5
V = range(n)
random.seed(int(12345))

all_possible = [frozenset(c) for k in [2,3,4] for c in itertools.combinations(V, k)]

def is_antichain(faces):
    for i in range(len(faces)):
        for j in range(len(faces)):
            if i != j and faces[i] <= faces[j]:
                return False
    return True

mismatches = []
tested_complexes = 0
tested_pairs = 0
q = 2

trials = 4000
for _ in range(trials):
    r = int(random.choice([2,3,3,4]))  # bias toward more faces
    combo = random.sample(all_possible, r)
    if not is_antichain(combo):
        continue
    covered = frozenset().union(*combo)
    if len(covered) < 5:
        continue  # skip if not all 5 vertices used (avoid degenerate free-vertex cases)
    maximal = [tuple(sorted(f)) for f in combo]
    Delta = complex_from_maximal(maximal)
    tested_complexes += 1
    for F in combo:
        if len(F) < 3:
            continue
        for e in itertools.combinations(sorted(F), 2):
            e = frozenset(e)
            tested_pairs += 1
            iso = is_isolating(Delta, F, e)
            ok, diff, chi_con = check_theorem(Delta, F, e, n, q)
            if iso != ok:
                mismatches.append((maximal, sorted(F), sorted(e), iso, ok))

print(f"complexes tested: {tested_complexes}, (F,e) pairs tested: {tested_pairs}")
print(f"mismatches: {len(mismatches)}")
for m in mismatches[:30]:
    print("  MISMATCH:", m)
