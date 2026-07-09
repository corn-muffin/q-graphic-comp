from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangements
import random

# Take a (0-indexed) simplicial complex delta and prime power q, and return the corresponding q-analogue arrangement
def q_graphic_arrangement(delta,q):
    n = len(delta.vertices())

    if delta.vertices() != tuple([0..n-1]):
        raise ValueError("Simplicial complex must be 0-indexed!")

    K = GF(q,'a')
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
    skeleton = delta.graph()
    clique_complex = skeleton.clique_complex()
    maybe = clique_complex.facets().list()
    #maybe = delta.graph().clique_complex().facets().list()
    facets.extend([f for f in maybe if random.random() <= p])
    return SimplicialComplex(facets)

# A: arrangement in GF(q)^A.ngens()
def sage_to_m2_arrangement(A,q):
    n= A.parent().ngens()

    # turn our hyperplanes into something m2 can read
    H = A.hyperplanes()
    m = len(H)
    M = "{"
    for j in [0..n-1]:
        row = "{"
        for i in [0..m-1]:
            #row += "a^("+ str(H[i][j]) +")"
            row += str(H[i][j])

            # we only want to separate inside a row
            if i<m-1:
                row += ", "
            else:
                row += "}"
        M += row

        # same exact idea
        if j<n-1:
            M+= ", "
        else:
            M+= "}"

    # turn it into m2 code
    m2_code = (
        f'needsPackage "HyperplaneArrangements";'
        +f' K = GF({q},Variable=> a);'
        +f' R = K[apply(0..{n-1},i->x_i)];'
        +f' A = arrangement(matrix{str(M)}, R)'
    )
    return m2_code

# take a hyperplane A over F_q and
# compute its projective dimension in m2
def compute_proj_dim(A,q):
    K= GF(q,'a')
    A.change_ring(K) # make sure the variable is in a!
    code_for_m2 = sage_to_m2_arrangement(A,q)
    code_for_m2 += "; pdim image der A" #find the pdim!
    return macaulay2.eval(code_for_m2) 

# take a simplex delta and compute the projective 
# dimension of S_delta^q in m2
# def compute_proj_dim(delta,q):
#     A = q_graphic_arrangement(delta,q)
#     return compute_proj_dim(A,q)
          

n = 5 # Number of vertices
d = 1 # Dimension of complex
p1 = 0.5 # Probability of including a given edge in 1-skeleton
p2 = 0.5 # Probability of including a face from the clique complex
s = 10 # Number of 1-skeleta to generate
t = 1 # Number of extensions of a given 1-skeleton
q = 2 # The finite field.
include_skeleta = 1

skeleta = sample_complexes(n,d,p1,s)

#print("working1")
num_cxs = s*(t+include_skeleta)
total_correct = 0
progress = 1
skeleta_tested = 0
counter_examples = []
proj_dims = []
for G in skeleta:
    #print("######################################### (" + str(progress) + "/" + str(num_cxs) + ")")

    #print("working2")
    complexes = complexes_on_skeleta(G,p2,t)
    if include_skeleta:
        complexes.append(G)
    #print("working3")
    #complexes = [extend_complex(G,p2)]
    w = G.graph().is_weakly_chordal()
    c = G.graph().is_chordal()

    for delta in complexes:
        print("######################################### (" + str(progress) + "/" + str(num_cxs) + ")")

        print(delta.facets())
        print("G weakly chordal: " + str(w))
        print("Delta is flag: " + str(delta.is_flag_complex()))
        A = q_graphic_arrangement(delta,q)
        pd = int(compute_proj_dim(A,q))
        proj_dims.append(pd)
        print("its projective dimension is: "+str(pd))
        if (pd<=1 and w) or (pd>1 and not w):
            print("Exactly what we expect!")
            total_correct+=1
        else:
            print("Something is wrong!")
            counter_examples.append(delta)

        progress += 1
        print("")

    #progress += 1
    print("")
