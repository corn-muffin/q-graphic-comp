load("qarr.sage")
import itertools

n = 5
maximal = [(0,1,2),(0,1,3),(0,1,4)]
Delta = complex_from_maximal(maximal)
F = frozenset({0,1,2})

Ddel = deletion(Delta, F)

subfaces = [frozenset(s) for k in [2,3] for s in itertools.combinations(sorted(F), k)]
print("subfaces of F (size>=2):", [sorted(s) for s in subfaces])

for q in [2,3,5]:
    chi_D = fast_characteristic_polynomial(Delta, n, q)
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
    diff = chi_D - chi_Ddel
    print(f"\n=== q={q} ===")
    print("diff =", diff.factor())
    for G in subfaces:
        Dcon, n_con = relabel_contiguous(contraction(Delta, G))
        chi_con = fast_characteristic_polynomial(Dcon, n_con, q)
        print(f"  G={sorted(G)}: chi(Delta/G) = {chi_con.factor()}   [divides diff? {chi_con.divides(diff)}]")
