load("qarr.sage")
import itertools, time

n = 5
V = range(n)
q = 2  # fast, just searching for ANY failure

def is_antichain(faces):
    faces = [frozenset(f) for f in faces]
    for i in range(len(faces)):
        for j in range(len(faces)):
            if i != j and faces[i] <= faces[j]:
                return False
    return True

candidates = []
sizes = [2, 3, 4]
all_possible = [frozenset(c) for k in sizes for c in itertools.combinations(V, k)]

t0 = time.time()
tested = 0
found = []
# try all pairs and triples of maximal faces forming antichains covering
# at least 4 of the 5 vertices (so the complex is "connected enough")
for r in [2, 3]:
    for combo in itertools.combinations(all_possible, r):
        if not is_antichain(combo):
            continue
        covered = frozenset().union(*combo)
        if len(covered) < 5:
            continue
        Delta = complex_from_maximal([tuple(sorted(f)) for f in combo])
        for F in combo:
            if len(F) < 3:
                continue  # edges always work by Nian's theorem; only test faces of dim>=2
            tested += 1
            Ddel = deletion(Delta, F)
            Dcon_raw = contraction(Delta, F)
            Dcon, n_con = relabel_contiguous(Dcon_raw)
            chi_D = fast_characteristic_polynomial(Delta, n, q)
            chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
            chi_Dcon = fast_characteristic_polynomial(Dcon, n_con, q)
            diff = chi_D - chi_Ddel
            if not chi_Dcon.divides(diff):
                found.append((tuple(sorted(map(sorted, combo))), tuple(sorted(F)), str(chi_D), str(chi_Ddel), str(chi_Dcon), str(diff/chi_Dcon)))
                print("FAILURE:", combo, "F=", sorted(F))
                print("   chi_D   =", chi_D.factor())
                print("   chi_Ddel=", chi_Ddel.factor())
                print("   chi_con =", chi_Dcon.factor())
                print("   diff/chi_con =", diff/chi_Dcon)
        if len(found) >= 8:
            break
    if len(found) >= 8:
        break

print(f"\ntested {tested} (complex, F) pairs in {time.time()-t0:.1f}s, found {len(found)} failures")
