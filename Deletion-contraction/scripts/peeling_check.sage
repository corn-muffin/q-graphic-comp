load("isolating_utils.sage")

n = 6
q = 2
maxA = [(0,1,2),(3,4,5),(0,3),(1,4),(2,5)]
Delta = complex_from_maximal(maxA)
F = frozenset({0,1,2})

# Step 1: Nian's PROVEN theorem on the maximal edge {0,3} (unconditionally valid)
e0 = frozenset({0,3})
Ddel_e0 = deletion(Delta, e0)
Dcon_e0, ncon_e0 = relabel_contiguous(contraction(Delta, e0))
chi_D = fast_characteristic_polynomial(Delta, n, q)
chi_Ddel_e0 = fast_characteristic_polynomial(Ddel_e0, n, q)
chi_Dcon_e0 = fast_characteristic_polynomial(Dcon_e0, ncon_e0, q)
nian_check = (chi_D == chi_Ddel_e0 - (q-1)*chi_Dcon_e0)
print("Step 1 (Nian's proven thm on bridge {0,3}) holds:", nian_check)

# Step 2: in Delta \ {0,3}, is F's edge {0,1} now isolating?
print("In Delta\\{0,3}: edge {0,1} of F isolating?", is_isolating(Ddel_e0, F, frozenset({0,1})))
print("In Delta\\{0,3}: edge {0,2} of F isolating?", is_isolating(Ddel_e0, F, frozenset({0,2})))

# Step 3: apply the (conjectured) isolating-edge theorem to Delta\{0,3} using e={0,1}
e = frozenset({0,1})
Ddel_e0_F = deletion(Ddel_e0, F)
Dcon_e, ncon_e = relabel_contiguous(contraction(Ddel_e0, e))
chi_Ddel_e0_F = fast_characteristic_polynomial(Ddel_e0_F, n, q)
chi_Dcon_e = fast_characteristic_polynomial(Dcon_e, ncon_e, q)
isolating_check = (chi_Ddel_e0 == chi_Ddel_e0_F - (q-1)**(len(F)-1)*chi_Dcon_e)
print("Step 2 (isolating thm applied inside Delta\\{0,3}) holds:", isolating_check)

# Combine: does chaining both steps reconstruct chi(Delta) exactly?
chained = chi_Ddel_e0_F - (q-1)**(len(F)-1)*chi_Dcon_e - (q-1)*chi_Dcon_e0
print("Chained formula reconstructs chi(Delta) exactly:", chained == chi_D)
