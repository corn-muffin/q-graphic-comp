from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangements
import random

# Take a (0-indexed) simplicial complex delta and prime power q, and return the corresponding q-analogue arrangement
def q_graphic_arrangement(delta,q):
    n = len(delta.vertices())

    if delta.vertices() != tuple([0..n-1]):
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

def complexes_on_skeleta(G,p,s):
    complexes = []

    while len(complexes) < s:
        delta = extend_complex(G,p)
        if delta not in complexes:
            complexes.append(delta)
    
    return complexes

def extend_complex(delta,p):
    vertices = delta.vertices()
    facets = delta.facets().list()
    maybe = delta.graph().clique_complex().facets().list()
    facets.extend([f for f in maybe if random.random() <= p])
    return SimplicialComplex(facets)

# A: arrangement in GF(q)^A.ngens()
def sage_to_m2_arrangement(A,q):
    #K= A.base_ring() # always the base field
    n= A.parent().ngens()

    # Extract polynomial strings
    H = A.hyperplanes()
    m = len(H)
    eqs = "{"
    for h in A.hyperplanes():
        row = "{"
        for i in [0..n-1]:
            eq += "a^("+ str(h[i]) + ")*x_" + str(i) +" + "
        eq += "a^("+str(h[i])+")"
        eqs.append(eq)

    # Generate Macaulay2 code blocks
    m2_code = (
        f'needsPackage "HyperplaneArrangements"\n'
        +f' K = GF({n},Variable=> a)\n'
        +f' R = K[apply(0..{n-1},i->x_i)]\n'
        +f' A = arrangement({str(eqs)}, K)'
    )
    return m2_code

# def compute_proj_dim(A,K)
# 
# def compute_proj_dim(delta,q):
#     A = q_graphic_arrangement(delta,q)
#     return compute_proj_dim(A,GF(q))
          

n = 4 # Number of vertices
d = 1 # Dimension of complex
p = 0.3 # Probability of including a face from the clique complex
s = 5 # Number of 1-skeleta to generate
q = 2

skeleta = sample_complexes(n,d,p,s)

progress = 1
for G in skeleta:
    print("######################################### (" + str(progress) + "/" + str(s) + ")")

    complexes = complexes_on_skeleta(G,p)
    print("G weakly chordal: " + str(G.graph().is_weakly_chordal()))

    for delta in complexes:
        print("######################################### (" + str(progress) + "/" + str(s) + ")")

        print(what

        progress += 1
        print("")

    progress += 1
    print("")
