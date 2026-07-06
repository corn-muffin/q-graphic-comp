load("qarr.sage")
n = 5
maximal = [(0,1,2),(0,1,3),(0,1,4),(2,3,4)]
Delta = complex_from_maximal(maximal)
F = frozenset({0,1,2})
Ddel = deletion(Delta, F)
for q in [2,3]:
    chi_D = fast_characteristic_polynomial(Delta, n, q)
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
    diff = chi_D - chi_Ddel
    e = frozenset({0,2})  # private (0 not in {2,3,4})
    Dcon, n_con = relabel_contiguous(contraction(Delta, e))
    chi_con = fast_characteristic_polynomial(Dcon, n_con, q)
    ok = (diff == -(q-1)**2 * chi_con)
    print(f"q={q}: private-edge thm holds for e={sorted(e)}? {ok}   (chi(Delta/e)={chi_con.factor()})")
