load("qarr.sage")
import itertools

n = 5
maximal = [(0,1,2,3),(0,1,2,4)]
Delta = complex_from_maximal(maximal)
F = frozenset({0,1,2,3})
Ddel = deletion(Delta, F)

edges = [frozenset(s) for s in itertools.combinations(sorted(F), 2)]

for q in [2,3]:
    chi_D = fast_characteristic_polynomial(Delta, n, q)
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
    diff = chi_D - chi_Ddel
    print(f"\n=== q={q} ===  diff = {diff.factor()}   (q-1)^(|F|-1) = {(q-1)**3}")
    for e in edges:
        Dcon, n_con = relabel_contiguous(contraction(Delta, e))
        chi_con = fast_characteristic_polynomial(Dcon, n_con, q)
        ratio_ok = (diff == -(q-1)**3 * chi_con)
        shared = "SHARED" if e <= frozenset({0,1,2}) else "private"
        print(f"  e={sorted(e)} ({shared}): chi(Delta/e) = {chi_con.factor()}   [diff == -(q-1)^3*chi_con ? {ratio_ok}]")
