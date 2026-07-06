load("qarr.sage")
n = 5
maximal = [(0,1,2,3),(0,1,2,4)]
Delta = complex_from_maximal(maximal)
F = frozenset({0,1,2,3})
Ddel = deletion(Delta, F)
private_ridges = [frozenset({0,1,3}), frozenset({0,2,3}), frozenset({1,2,3})]  # size-3, private
private_edges = [frozenset({0,3}), frozenset({1,3}), frozenset({2,3})]

for q in [2,3]:
    chi_D = fast_characteristic_polynomial(Delta, n, q)
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
    diff = chi_D - chi_Ddel
    print(f"\n=== q={q} ===  diff = {diff.factor()}")
    for G in private_ridges:
        Dcon, n_con = relabel_contiguous(contraction(Delta, G))
        chi_con = fast_characteristic_polynomial(Dcon, n_con, q)
        for coeff, label in [(-(q-1)**3, "(q-1)^3"), (-(q-1)**2, "(q-1)^2"), (-(q-1), "(q-1)^1")]:
            if diff == coeff*chi_con:
                print(f"  ridge G={sorted(G)}: chi(Delta/G)={chi_con.factor()}  MATCHES with coeff {label}")
                break
        else:
            print(f"  ridge G={sorted(G)}: chi(Delta/G)={chi_con.factor()}  divides diff? {chi_con.divides(diff)}  diff/chi_con={diff/chi_con if chi_con.divides(diff) else 'N/A'}")
