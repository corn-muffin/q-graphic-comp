load("qarr.sage")
import time

cases = {
    "full_simplex_4":      dict(n=4, maximal=[(0,1,2,3)], F=(0,1,2,3)),
    "boundary_simplex_4":  dict(n=4, maximal=[(0,1,2),(0,1,3),(0,2,3),(1,2,3)], F=(0,1,2)),
    "book_2pages":         dict(n=4, maximal=[(0,1,2),(0,1,3)], F=(0,1,2)),
    "book_3pages":         dict(n=5, maximal=[(0,1,2),(0,1,3),(0,1,4)], F=(0,1,2)),
    "book_3pages_plus_opposite_triangle": dict(n=5, maximal=[(0,1,2),(0,1,3),(0,1,4),(2,3,4)], F=(0,1,2)),
    "nocone_C":            dict(n=5, maximal=[(0,1,2),(2,3,4),(0,4)], F=(0,1,2)),
}

for name, spec in cases.items():
    n = spec["n"]
    Delta = complex_from_maximal(spec["maximal"])
    F = frozenset(spec["F"])
    Ddel = deletion(Delta, F)
    Dcon, n_con = relabel_contiguous(contraction(Delta, F))
    q = 5
    t0 = time.time()
    chi_D = fast_characteristic_polynomial(Delta, n, q)
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
    chi_Dcon = fast_characteristic_polynomial(Dcon, n_con, q)
    dt = time.time() - t0
    diff = chi_D - chi_Ddel
    div = chi_Dcon.divides(diff)
    print(f"{name}  F={sorted(F)}  (q=5, {dt:.1f}s)")
    print(f"   chi_D    = {chi_D.factor()}")
    print(f"   chi_Ddel = {chi_Ddel.factor()}")
    print(f"   chi_con  = {chi_Dcon.factor()}")
    print(f"   divisible={div}  quotient/frac = {(diff/chi_Dcon)}")
    print()
