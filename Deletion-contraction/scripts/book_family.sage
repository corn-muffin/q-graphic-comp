load("qarr.sage")

# "Book" of k triangles all sharing the edge {0,1}: pages {0,1,2},{0,1,3},...,{0,1,k+1}
# Uses n = k+2 vertices.  Test deletion-contraction on F = a single page {0,1,2}.
for k in [1, 2, 3, 4]:
    n = k + 2
    pages = [(0, 1, v) for v in range(2, k + 2)]
    Delta = complex_from_maximal(pages)
    F = frozenset({0, 1, 2})
    Ddel = deletion(Delta, F)
    Dcon, n_con = relabel_contiguous(contraction(Delta, F))
    q = 2
    chi_D = fast_characteristic_polynomial(Delta, n, q)
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
    chi_Dcon = fast_characteristic_polynomial(Dcon, n_con, q)
    diff = chi_D - chi_Ddel
    div = chi_Dcon.divides(diff)
    print(f"k={k} pages (n={n}): chi_D={chi_D.factor()}")
    print(f"           chi_Ddel={chi_Ddel.factor()}")
    print(f"           chi_con ={chi_Dcon.factor()}")
    print(f"           divisible={div}  diff/chi_con = {(diff/chi_Dcon).factor() if not div else diff/chi_Dcon}")
    print()
