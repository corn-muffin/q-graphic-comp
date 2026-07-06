load("qarr.sage")

n = 7
# book_3pages plus an extra triangle {1,5,6} attached only to spine vertex 1
# (definitively breaks 0<->1 swap symmetry)
maximal = [(0,1,2),(0,1,3),(0,1,4),(1,5,6)]
Delta = complex_from_maximal(maximal)
F = frozenset({0,1,2})
Ddel = deletion(Delta, F)
private_edges = [frozenset({0,2}), frozenset({1,2})]

for q in [2,3]:
    chi_D = fast_characteristic_polynomial(Delta, n, q)
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
    diff = chi_D - chi_Ddel
    print(f"\n=== q={q} ===  diff = {diff.factor()}")
    for e in private_edges:
        Dcon, n_con = relabel_contiguous(contraction(Delta, e))
        chi_con = fast_characteristic_polynomial(Dcon, n_con, q)
        ok = (diff == -(q-1)**2 * chi_con)
        print(f"  e={sorted(e)}: chi(Delta/e) = {chi_con.factor()}   [thm holds? {ok}]")
