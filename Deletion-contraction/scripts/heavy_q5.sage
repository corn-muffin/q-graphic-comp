load("qarr.sage")
import time, json

results = {}
targets = [
    ("full_simplex_5", 5, [(0,1,2,3,4)], (0,1,2,3,4)),
    ("boundary_simplex_5", 5, [(0,1,2,3),(0,1,2,4),(0,1,3,4),(0,2,3,4),(1,2,3,4)], (0,1,2,3)),
]
q = 5
for name, n, maximal, F in targets:
    Delta = complex_from_maximal(maximal)
    F = frozenset(F)
    Ddel = deletion(Delta, F)
    Dcon, n_con = relabel_contiguous(contraction(Delta, F))
    print(f"=== {name} (q=5) ===", flush=True)
    t0 = time.time()
    chi_D = fast_characteristic_polynomial(Delta, n, q, verbose=True)
    print(f"  chi_D done in {time.time()-t0:.1f}s: {chi_D.factor()}", flush=True)
    t0 = time.time()
    chi_Ddel = fast_characteristic_polynomial(Ddel, n, q, verbose=True)
    print(f"  chi_Ddel done in {time.time()-t0:.1f}s: {chi_Ddel.factor()}", flush=True)
    t0 = time.time()
    chi_Dcon = fast_characteristic_polynomial(Dcon, n_con, q, verbose=True)
    print(f"  chi_Dcon done in {time.time()-t0:.1f}s: {chi_Dcon.factor()}", flush=True)
    diff = chi_D - chi_Ddel
    div = chi_Dcon.divides(diff)
    print(f"  divisible={div}  quotient/frac={diff/chi_Dcon}", flush=True)
    results[name] = dict(chi_D=str(chi_D), chi_Ddel=str(chi_Ddel), chi_Dcon=str(chi_Dcon),
                          divisible=bool(div), quotient=str(diff/chi_Dcon))
    with open("../data/heavy_q5_results.json", "w") as f:
        json.dump(results, f, indent=2)
print("ALL DONE")
