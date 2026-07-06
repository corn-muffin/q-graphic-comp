load("qarr.sage")
n = 6
q = 2
F = frozenset({0,1,2})

# Prism sides open, but WITHOUT the {3,4,5} triangle (just 3 disjoint pendant edges)
maxD = [(0,1,2),(0,3),(1,4),(2,5)]
DeltaD = complex_from_maximal(maxD)
DdelD = deletion(DeltaD, F)
chi_D = fast_characteristic_polynomial(DeltaD, n, q)
chi_Ddel = fast_characteristic_polynomial(DdelD, n, q)
diff = chi_D - chi_Ddel
print("3 independent pendants (no {3,4,5} face), diff =", diff.factor())
for e in [frozenset({0,1}), frozenset({0,2})]:
    Dcon, nc = relabel_contiguous(contraction(DeltaD, e))
    chi_con = fast_characteristic_polynomial(Dcon, nc, q)
    ok = (diff == -(q-1)**2*chi_con)
    print(f"  e={sorted(e)}: chi(Delta/e)={chi_con.factor()}  thm_holds={ok}")

# Prism with {3,4,5} triangle but only ONE connecting edge {0,3} (not all 3)
maxE = [(0,1,2),(3,4,5),(0,3)]
DeltaE = complex_from_maximal(maxE)
DdelE = deletion(DeltaE, F)
chi_D2 = fast_characteristic_polynomial(DeltaE, n, q)
chi_Ddel2 = fast_characteristic_polynomial(DdelE, n, q)
diff2 = chi_D2 - chi_Ddel2
print("\n{3,4,5} triangle + only ONE connecting edge {0,3}, diff =", diff2.factor())
for e in [frozenset({0,1}), frozenset({1,2})]:
    Dcon, nc = relabel_contiguous(contraction(DeltaE, e))
    chi_con = fast_characteristic_polynomial(Dcon, nc, q)
    ok = (diff2 == -(q-1)**2*chi_con)
    print(f"  e={sorted(e)}: chi(Delta/e)={chi_con.factor()}  thm_holds={ok}")
