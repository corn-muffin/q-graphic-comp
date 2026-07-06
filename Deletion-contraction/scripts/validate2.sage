load("qarr.sage")
import time

Delta53 = skeleton(full_simplex(5), 3)

expected = {
    2: [1,2,4,8,10],
    3: [1,3,9,25,27],
    5: [1,5,25],  # plus irreducible quadratic t^2-174t+7661
}

print("=== Fast method validation on Delta_{5,3} ===")
for q in [2, 3, 5]:
    t0 = time.time()
    chi = fast_characteristic_polynomial(Delta53, 5, q, verbose=True)
    dt = time.time() - t0
    print(f"q={q}: chi = {chi.factor()}")
    print(f"   time = {dt:.2f}s")
