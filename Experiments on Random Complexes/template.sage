from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangements
import random

# Take a (0-indexed) simplicial complex delta and prime power q, and return the corresponding q-analogue arrangement
def q_graphic_arrangement(delta,q):
    n = len(delta.vertices())

    if delta.vertices() != [0..n-1]:
        raise ValueError("Simplicial complex must be 0-indexed!")

    K = GF(q)
    H = HyperplaneArrangements(K,tuple('x'+str(i) for i in [0..n-1]))
    x = H.gens()

    facets = delta.maximal_faces()

    A = H()
    for face in facets:
        m = dim(face) + 1 # Number of vertices of face
        V = VectorSpace(K,m) # m-tuples of coefficients for hyperplanes
        
        for v in V:
            if v != V.zero_vector():
                j = 0
                a = [] # Coefficients

                # Insert coefficients from vector to entries of a indexed by face
                for i in [0..n-1]:
                    if i in face:
                        a.append(v[j])
                        j += 1
                    else:
                        a.append(0)
                
                A = A.add_hyperplane(sum(a[i]*x[i] for i in [0..n-1]))
    return A


# Compute characteristic polynomial of q-graphic arrangement associated to a simplicial complex delta for fixed q
def q_graphic_char_polynomial(delta,q):
    S = q_graphic_arrangement(delta,q)

    # characteristic_polynomial() doesn't work when q isn't prime
    try:
        return S.characteristic_polynomial()
    except Exception:
        return S._slow_characteristic_polynomial()


# Generate s simplicial complexes on n vertices of dimension d
def sample_complexes(n,d,p,s):
    complexes = []

    while len(complexes) < s:
        delta = simplicial_complexes.RandomComplex(n,d,p)
        if delta not in complexes:
            complexes.append(delta)
    
    return complexes


# Compute contraction of delta by F (courtesy of Claude; should be double-checked)
def contract_simplicial_complex(delta, F):
    F_verts = frozenset(F)
    facets = [frozenset(f) for f in delta.facets()]

    # Choose a representative vertex for the identified point
    v_star = min(F_verts)

    def phi(v):
        return v_star if v in F_verts else v

    # Push every facet forward under phi.
    new_faces = {frozenset(phi(v) for v in facet) for facet in facets}

    maximal_faces = [f for f in new_faces
                      if not any(f < g for g in new_faces if f != g)]

    # Reindex vertices to 0, ..., n-1, ordered by their original labels
    all_verts = sorted(set().union(*maximal_faces))
    relabel = {v: i for i, v in enumerate(all_verts)}

    reindexed_faces = [[relabel[v] for v in f] for f in maximal_faces]

    return SimplicialComplex(reindexed_faces)


n = 6 # Number of vertices
d = 4 # Dimension of complex
p = 0.5 # Probability of including a d-face
s = 10 # Number of simplicial complexes to generate

complexes = sample_complexes(n,d,p,s)

progress = 1
for delta in complexes:
    print("######################################### (" + str(progress) + "/" + str(s) + ")")

    # [Procedure here]

    progress += 1
    print("")