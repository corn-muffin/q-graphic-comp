load("qarr.sage")
n = 6
maxA = [(0,1,2),(3,4,5),(0,3),(1,4),(2,5)]
DeltaA = complex_from_maximal(maxA)
F = frozenset({0,1,2})
Ddel = deletion(DeltaA, F)
q = 2
chi_D = fast_characteristic_polynomial(DeltaA, n, q)
chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
diff = chi_D - chi_Ddel
print("diff =", diff.factor())
for e in [frozenset({0,1}), frozenset({0,2}), frozenset({1,2})]:
    Dcon, nc = relabel_contiguous(contraction(DeltaA, e))
    chi_con = fast_characteristic_polynomial(Dcon, nc, q)
    ok = (diff == -(q-1)**2*chi_con)
    print(f"  e={sorted(e)}: chi(Delta/e)={chi_con.factor()}  thm_holds={ok}")

# Now compare: book_3pages but with vertex 2 ALSO given a pendant (breaking its isolation)
maxC = [(0,1,2),(0,1,3),(0,1,4),(2,5)]
DeltaC = complex_from_maximal(maxC)
n2 = 6
Fc = frozenset({0,1,2})
DdelC = deletion(DeltaC, Fc)
chi_Dc = fast_characteristic_polynomial(DeltaC, n2, q)
chi_DdelC = fast_characteristic_polynomial(DdelC, n2, q)
diffC = chi_Dc - chi_DdelC
print("\nbook_3pages + pendant at vertex 2 (breaks its isolation): diff =", diffC.factor())
for e in [frozenset({0,2}), frozenset({1,2})]:
    Dcon, nc = relabel_contiguous(contraction(DeltaC, e))
    chi_con = fast_characteristic_polynomial(Dcon, nc, q)
    ok = (diffC == -(q-1)**2*chi_con)
    print(f"  e={sorted(e)}: chi(Delta/e)={chi_con.factor()}  thm_holds={ok}")
