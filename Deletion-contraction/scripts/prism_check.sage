load("qarr.sage")

n = 6
# Version A: two triangles {0,1,2},{3,4,5} + matching edges 0-3,1-4,2-5 (sides NOT filled)
maxA = [(0,1,2),(3,4,5),(0,3),(1,4),(2,5)]
DeltaA = complex_from_maximal(maxA)

# Version B: fully triangulated closed prism surface (sides filled with 2 triangles each)
maxB = [(0,1,2),(3,4,5),(0,1,4),(0,4,3),(1,2,5),(1,5,4),(2,0,3),(2,3,5)]
DeltaB = complex_from_maximal(maxB)

F = frozenset({0,1,2})

for name, Delta in [("A: prism, sides open", DeltaA), ("B: prism, sides triangulated (pseudomanifold)", DeltaB)]:
    print(f"\n=== {name} ===")
    edges = [frozenset(s) for s in [(0,1),(0,2),(1,2)]]
    for e in edges:
        others = [G for G in Delta if len(G)>=2 and e <= G and not (G <= F)]
        print(f"  edge {sorted(e)} private? {len(others)==0}  (other faces containing it beyond F: {[sorted(g) for g in others]})")
    Ddel = deletion(Delta, F)
    for q in [2,3]:
        chi_D = fast_characteristic_polynomial(Delta, n, q)
        chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
        diff = chi_D - chi_Ddel
        Dconf, ncf = relabel_contiguous(contraction(Delta, F))
        chi_conf = fast_characteristic_polynomial(Dconf, ncf, q)
        divF = chi_conf.divides(diff)
        print(f"  q={q}: [Delta/F direct] divides={divF}  p_F={(diff/chi_conf) if divF else 'N/A'}")
        e = frozenset({0,1})
        Dcone, nce = relabel_contiguous(contraction(Delta, e))
        chi_cone = fast_characteristic_polynomial(Dcone, nce, q)
        ok = (diff == -(q-1)**2 * chi_cone)
        print(f"  q={q}: [private-edge thm via {{0,1}}] holds={ok}  chi(Delta/e)={chi_cone.factor()}")
