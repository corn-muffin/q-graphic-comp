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

danger = []       # iso=True but ok=False -- REAL counterexample to sufficiency
overperform = 0    # iso=False but ok=True -- isolating just not necessary
tested_pairs = 0
q = 2
trials = 4000
for _ in range(trials):
    r = int(random.choice([2,3,3,4]))
    combo = random.sample(all_possible, r)
    if not is_antichain(combo):
        continue
    covered = frozenset().union(*combo)
    if len(covered) < 5:
        continue
    maximal = [tuple(sorted(f)) for f in combo]
    Delta = complex_from_maximal(maximal)
    for F in combo:
        if len(F) < 3:
            continue
        for e in itertools.combinations(sorted(F), 2):
            e = frozenset(e)
            tested_pairs += 1
            iso = is_isolating(Delta, F, e)
            ok, diff, chi_con = check_theorem(Delta, F, e, n, q)
            if iso and not ok:
                danger.append((maximal, sorted(F), sorted(e)))
            elif not iso and ok:
                overperform += 1

print(f"(F,e) pairs tested: {tested_pairs}")
print(f"DANGEROUS mismatches (iso=True but theorem fails): {len(danger)}")
for d in danger[:20]:
    print("  ", d)
print(f"Overperformance (iso=False but theorem still holds anyway): {overperform}")
