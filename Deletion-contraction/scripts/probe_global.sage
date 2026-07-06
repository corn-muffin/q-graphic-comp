load("qarr.sage")

n = 5
q = 2
base_pages = [(0,1,2),(0,1,3),(0,1,4)]  # the failing k=3 book
F = frozenset({0,1,2})

def test(name, extra_pages):
    maximal = base_pages + extra_pages
    Delta = complex_from_maximal(maximal)
    Ddel = deletion(Delta, F)
    Dcon, n_con = relabel_contiguous(contraction(Delta, F))
    chi_D = fast_characteristic_polynomial(Delta, n, q)
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
    chi_Dcon = fast_characteristic_polynomial(Dcon, n_con, q)
    diff = chi_D - chi_Ddel
    div = chi_Dcon.divides(diff)
    print(f"{name}: maximal={maximal}")
    print(f"   divisible={div}   diff/chi_con = {diff/chi_Dcon}")
    print()

test("baseline (fails)", [])
test("+ triangle {2,3,4} disjoint from edge {0,1}", [(2,3,4)])
test("+ triangle {0,2,3} sharing vertex 0 only", [(0,2,3)])
test("+ triangle {1,2,3} sharing vertex 1 only", [(1,2,3)])
test("+ edge {2,3} only (no new triangle)", [(2,3)])
test("+ all of {0,2,3},{0,2,4},{0,3,4}", [(0,2,3),(0,2,4),(0,3,4)])
test("+ all of {1,2,3},{1,2,4},{1,3,4}", [(1,2,3),(1,2,4),(1,3,4)])
test("+ {2,3,4} and {0,2,3}", [(2,3,4),(0,2,3)])
test("full 2-skeleton (all 10 triangles)", [(0,2,3),(0,2,4),(0,3,4),(1,2,3),(1,2,4),(1,3,4),(2,3,4)])
