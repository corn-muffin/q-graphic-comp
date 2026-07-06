load("qarr.sage")

# Delta_{5,3} = 2-skeleton of the 4-simplex on 5 vertices (faces of size <= 3).
Delta53 = skeleton(full_simplex(5), 3)

print("=== Validating against Nian's thesis Example 3.5(2), Delta_{5,3} ===")
for q in [2, 3, 4, 5]:
    chi = characteristic_polynomial(Delta53, 5, q, verbose=True)
    print(f"q={q}: chi = {chi.factor()}")
    print(f"       = {chi}")
