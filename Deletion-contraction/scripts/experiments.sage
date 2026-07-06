load("qarr.sage")
import time, json, sys

# ---------------------------------------------------------------------------
# Test complexes.  All on vertex set {0,...,n-1}.  Given by their maximal
# faces (tuples of vertex indices); complex_from_maximal takes the downward
# closure.
# ---------------------------------------------------------------------------

COMPLEXES = {
    # --- baseline family from Nian's thesis (sanity checks / new data) ---
    "full_simplex_4":      dict(n=4, maximal=[(0,1,2,3)]),
    "boundary_simplex_4":  dict(n=4, maximal=[(0,1,2),(0,1,3),(0,2,3),(1,2,3)]),
    "full_simplex_5":      dict(n=5, maximal=[(0,1,2,3,4)]),
    "boundary_simplex_5":  dict(n=5, maximal=[(0,1,2,3),(0,1,2,4),(0,1,3,4),(0,2,3,4),(1,2,3,4)]),
    "skeleton_5_3":        dict(n=5, maximal=[c for c in __import__("itertools").combinations(range(5),3)]),

    # --- cone-point test cases (each has a maximal 2-face / triangle F) ---
    # Case A: cone point 0 is in EVERY maximal face.
    #   maximal faces: {0,1,2,3} (tetrahedron) and {0,1,4} (triangle)
    "cone_A":              dict(n=5, maximal=[(0,1,2,3),(0,1,4)]),
    # Case B: same shape, cone point 0, but F is the "small" triangle sharing
    #   an edge {0,1} with the big facet (structurally like case A) -- kept
    #   as a duplicate-style sanity check with different geometry:
    #   maximal faces: {0,1,2,3} and {0,3,4}
    "cone_B":              dict(n=5, maximal=[(0,1,2,3),(0,3,4)]),
    # Case C: NO cone point (no vertex lies in all maximal faces).
    #   maximal faces: {0,1,2} (triangle), {2,3,4} (triangle), {0,4} (edge)
    "nocone_C":            dict(n=5, maximal=[(0,1,2),(2,3,4),(0,4)]),
    # Case D: NO cone point, all maximal faces triangles (pure 2-dim, no
    #   apex): {0,1,2}, {2,3,4}, {4,0,1}  -- a "triangulated cycle"
    "nocone_D":            dict(n=5, maximal=[(0,1,2),(2,3,4),(4,0,1)]),
    # Case E: cone point 0, base has two disjoint maximal edges (no shared
    #   vertex in the base): maximal faces {0,1,2} and {0,3,4}
    "cone_E":              dict(n=5, maximal=[(0,1,2),(0,3,4)]),

    # --- the "book" family: k triangles glued along a common edge {0,1} ---
    "book_2pages":         dict(n=4, maximal=[(0,1,2),(0,1,3)]),
    "book_3pages":         dict(n=5, maximal=[(0,1,2),(0,1,3),(0,1,4)]),  # MINIMAL COUNTEREXAMPLE
    "book_3pages_plus_opposite_triangle": dict(n=5, maximal=[(0,1,2),(0,1,3),(0,1,4),(2,3,4)]),
    "book_3pages_plus_cone0": dict(n=5, maximal=[(0,1,2),(0,1,3),(0,1,4),(0,2,3),(0,2,4),(0,3,4)]),
}


def analyze_face(Delta, F, n, qs, log):
    """For maximal face F of Delta, compute chi(Delta), chi(Delta\\F),
    chi(Delta/F) at each q in qs, and check polynomial divisibility of
    the difference by chi(Delta/F)."""
    F = frozenset(F)
    Ddel = deletion(Delta, F)
    Dcon_raw = contraction(Delta, F)
    Dcon, n_con = relabel_contiguous(Dcon_raw)

    result = {"F": sorted(F), "n_contracted": n_con, "per_q": {}}
    for q in qs:
        t0 = time.time()
        chi_D = fast_characteristic_polynomial(Delta, n, q)
        chi_Ddel = fast_characteristic_polynomial(Ddel, n, q)
        chi_Dcon = fast_characteristic_polynomial(Dcon, n_con, q)
        dt = time.time() - t0
        diff = chi_D - chi_Ddel
        R = diff.parent()
        divisible = chi_Dcon.divides(diff)
        entry = {
            "chi_Delta": str(chi_D),
            "chi_Delta_minus_F": str(chi_Ddel),
            "chi_Delta_over_F": str(chi_Dcon),
            "diff": str(diff),
            "divisible": bool(divisible),
            "time": round(dt, 2),
        }
        if divisible:
            quot = R(diff / chi_Dcon)
            entry["quotient_p_F"] = str(quot)
        else:
            # report as a reduced rational function for inspection
            frac = diff / chi_Dcon
            entry["quotient_p_F"] = None
            entry["diff_over_chi_con_as_fraction"] = str(frac)
        result["per_q"][q] = entry
        log(f"    q={q}: divisible={divisible}  "
            f"{'p_F='+entry['quotient_p_F'] if divisible else 'NOT POLY: diff/chi_con='+entry['diff_over_chi_con_as_fraction']}  ({dt:.2f}s)")
    return result


def main():
    qs = [2, 3]
    if len(sys.argv) > 1 and sys.argv[1] == "with5":
        qs = [2, 3, 5]

    all_results = {}
    for name, spec in COMPLEXES.items():
        n = spec["n"]
        Delta = complex_from_maximal(spec["maximal"])
        maxfaces = sorted(spec["maximal"], key=len, reverse=True)
        print(f"\n=== {name}  (n={n}, maximal faces={spec['maximal']}) ===")
        all_results[name] = {"n": n, "maximal": spec["maximal"], "faces": []}
        for F in maxfaces:
            if len(F) < 2:
                continue  # skip trivial vertex faces
            print(f"  -- testing deletion-contraction on maximal face F={sorted(F)} (|F|={len(F)}) --")
            res = analyze_face(Delta, F, n, qs, log=print)
            all_results[name]["faces"].append(res)

    def sagefix(o):
        if isinstance(o, dict):
            return {sagefix(k): sagefix(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [sagefix(x) for x in o]
        try:
            return int(o)
        except (TypeError, ValueError):
            return o

    with open("../data/experiment_results.json", "w") as f:
        json.dump(sagefix(all_results), f, indent=2)
    print("\nSaved to ../data/experiment_results.json")


main()
