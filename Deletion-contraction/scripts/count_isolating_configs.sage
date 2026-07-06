load("isolating_utils.sage")
import itertools

n = 5
V = range(n)

P = Posets.BooleanLattice(5)
# elements of P are integers 0..31, bit i set means vertex i is in the subset
def bits_to_frozenset(x):
    return frozenset(i for i in range(5) if (x >> i) & 1)

total_antichains = 0
total_complexes_5vertex = 0
total_connected = 0
count_any_F = 0       # exists maximal F (|F|>=2) with an isolating private edge
count_F_ge3 = 0        # exists maximal F with |F|>=3 with an isolating private edge

examples_ge3 = []

for A in P.antichains():
    total_antichains += 1
    faces = [bits_to_frozenset(x) for x in A if x != 0]   # drop the empty set
    if not faces:
        continue
    covered = frozenset().union(*faces)
    if len(covered) < 5:
        continue
    total_complexes_5vertex += 1
    maximal = [tuple(sorted(f)) for f in faces]
    Delta = complex_from_maximal(maximal)
    if not is_connected(Delta):
        continue
    total_connected += 1

    found_any = False
    found_ge3 = False
    for F in faces:
        if len(F) < 2:
            continue
        for e in itertools.combinations(sorted(F), 2):
            e = frozenset(e)
            if is_private(Delta, F, e) and is_isolating(Delta, F, e):
                found_any = True
                if len(F) >= 3:
                    found_ge3 = True
        if found_ge3:
            break

    if found_any:
        count_any_F += 1
    if found_ge3:
        count_F_ge3 += 1
        if len(examples_ge3) < 5:
            examples_ge3.append(maximal)

print("total antichains (Dedekind check):", total_antichains)
print("complexes spanning all 5 vertices:", total_complexes_5vertex)
print("...of those, connected:", total_connected)
print("connected & has SOME maximal F (any size>=2) with isolating private edge:", count_any_F)
print("connected & has SOME maximal F with |F|>=3 and an isolating private edge:", count_F_ge3)
print()
print("sample |F|>=3 examples:")
for ex in examples_ge3:
    print("  ", ex)
