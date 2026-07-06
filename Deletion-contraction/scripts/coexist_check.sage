load("qarr.sage")

n = 5
maximal = [(0,1,2,3),(0,1,4)]
Delta = complex_from_maximal(maximal)
F = frozenset({0,1,4})
Ddel = deletion(Delta, F)
private_edges = [frozenset({0,4}), frozenset({1,4})]  # {0,1} is shared with {0,1,2,3}

for q in [2,3]:
    chi_D = fast_characteristic_polynomial(Delta, n, q)
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
    diff = chi_D - chi_Ddel
    print(f"\n=== q={q} ===  diff = {diff.factor()}")
    Dcon_full, ncf = relabel_contiguous(contraction(Delta, F))
    chi_conF = fast_characteristic_polynomial(Dcon_full, ncf, q)
    print(f"  [old, full F] chi(Delta/F) = {chi_conF.factor()}, divides? {chi_conF.divides(diff)}, p_F={diff/chi_conF if chi_conF.divides(diff) else 'N/A'}")
    for e in private_edges:
        Dcon, n_con = relabel_contiguous(contraction(Delta, e))
        chi_con = fast_characteristic_polynomial(Dcon, n_con, q)
        ok = (diff == -(q-1)**2 * chi_con)
        print(f"  [new, private edge {sorted(e)}] chi(Delta/e) = {chi_con.factor()}   [thm holds? {ok}]")
