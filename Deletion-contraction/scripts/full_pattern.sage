load("qarr.sage")
n = 5
maximal = [(0,1,2,3),(0,1,2,4)]
Delta = complex_from_maximal(maximal)
F = frozenset({0,1,2,3})
Ddel = deletion(Delta, F)

for q in [2,3,5]:
    chi_D = fast_characteristic_polynomial(Delta, n, q)
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
    diff = chi_D - chi_Ddel
    Dcon, n_con = relabel_contiguous(contraction(Delta, F))
    chi_con = fast_characteristic_polynomial(Dcon, n_con, q)
    print(f"q={q}: chi(Delta/F)={chi_con.factor()}  divides={chi_con.divides(diff)}  quotient={(diff/chi_con).factor() if chi_con.divides(diff) else 'N/A'}")
