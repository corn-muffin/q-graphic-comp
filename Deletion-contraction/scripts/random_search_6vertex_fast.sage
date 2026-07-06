load("isolating_utils.sage")
import itertools, random

n = 6
V = range(n)
random.seed(int(20260705))
all_possible = [frozenset(c) for k in [2,3,4,5,6] for c in itertools.combinations(V, k)]

def is_antichain(faces):
    for i in range(len(faces)):
        for j in range(len(faces)):
            if i != j and faces[i] <= faces[j]:
                return False
    return True

danger = []
overperform = 0
matched_iso = 0
tested_pairs = 0
connected_complexes = 0
disconnected_skipped = 0
complexes_with_isolating_edge = 0
q = 2
trials = 1500

for trial_i in range(trials):
    r = int(random.choice([2,2,3,3,4]))
    combo = random.sample(all_possible, r)
    if not is_antichain(combo):
        continue
    covered = frozenset().union(*combo)
    if len(covered) < n:
        continue
    maximal = [tuple(sorted(f)) for f in combo]
    Delta = complex_from_maximal(maximal)
    if not is_connected(Delta):
        disconnected_skipped += 1
        continue
    connected_complexes += 1

    chi_D = fast_characteristic_polynomial(Delta, n, q)  # once per complex

    has_iso_here = False
    for F in combo:
        if len(F) < 3:
            continue
        F = frozenset(F)
        Ddel = deletion(Delta, F)
        chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)  # once per (Delta,F)
        diff = chi_D - chi_Ddel
        for e in itertools.combinations(sorted(F), 2):
            e = frozenset(e)
            tested_pairs += 1
            iso = is_isolating(Delta, F, e)
            Dcon, nc = relabel_contiguous(contraction(Delta, e))
            chi_con = fast_characteristic_polynomial(Dcon, nc, q)
            ok = (diff == -(q-1)**(len(F)-1) * chi_con)
            if iso:
                has_iso_here = True
                if ok:
                    matched_iso += 1
                else:
                    danger.append((maximal, sorted(F), sorted(e)))
            else:
                if ok:
                    overperform += 1
    if has_iso_here:
        complexes_with_isolating_edge += 1
    if (trial_i+1) % 200 == 0:
        print(f"...progress: {trial_i+1}/{trials} trials, {connected_complexes} connected so far", flush=True)

print(f"trials: {trials}")
print(f"connected complexes tested: {connected_complexes}")
print(f"disconnected skipped: {disconnected_skipped}")
print(f"complexes with >=1 isolating edge on some F(|F|>=3): {complexes_with_isolating_edge}")
print(f"(F,e) pairs tested: {tested_pairs}")
print(f"isolating & theorem holds (expected): {matched_iso}")
print(f"DANGEROUS mismatches (isolating but theorem FAILS): {len(danger)}")
for d in danger[:20]:
    print("  ", d)
print(f"overperformance (non-isolating, theorem holds anyway): {overperform}")
