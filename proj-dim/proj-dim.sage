from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangements
import random

# Take a (0-indexed) simplicial complex delta and prime power q, and return the corresponding q-analogue arrangement
def q_graphic_arrangement(delta,q):
    n = len(delta.vertices())

    if delta.vertices() != tuple([0..n-1]):
        raise ValueError("Simplicial complex must be 0-indexed!")

    K = QQ[
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

# Take a (0-indexed) simplicial complex delta and return the corresponding graphic arrangement
def graphic_arrangement(delta):
    n = len(delta.vertices())

    if delta.vertices() != tuple([0..n-1]):
        raise ValueError("Simplicial complex must be 0-indexed!")

    K = QQ
    H = HyperplaneArrangements(K,tuple('x'+str(i) for i in [0..n-1]))
    x = H.gens()

    facets = delta.maximal_faces()

    A = H()
    for face in facets:
        m = dim(face) + 1 # Number of vertices of face
        V = VectorSpace(K,m) # m-tuples of coefficients for hyperplanes

        A.add_hyperplane()

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
          

n = 8 # Number of vertices
d = 1 # Dimension of complex
p1 = 0.4 # Probability of including a given edge in 1-skeleton
p2 = 0.5 # Probability of including a face from the clique complex
s = 10 # Number of 1-skeleta to generate
t = 0 # Number of extensions of a given 1-skeleton
q = 2 # The finite field.
include_skeleta = 1
include_clique = 1

skeleta = sample_complexes(n,d,p1,s)

#print("working1")
num_cxs = s*(t+include_skeleta + include_clique)
really_bad = 0
progress = 1
skeleta_tested = 0
tested_complexes = []
counter_examples = []
proj_dims = []
for G in skeleta:
    #print("######################################### (" + str(progress) + "/" + str(num_cxs) + ")")

    #print("working2")
    complexes = complexes_on_skeleta(G,p2,t)
    if include_skeleta:
        complexes.append(G)
    if include_clique:
        complexes.append(G.graph().clique_complex())
    #print("working3")
    #complexes = [extend_complex(G,p2)]
    w = G.graph().is_weakly_chordal()
    c = G.graph().is_chordal()

    for delta in complexes:
        print("######################################### (" + str(progress) + "/" + str(num_cxs) + ")")

        print(delta.facets())
        print("G weakly chordal: " + str(w))
        print("G chordal: " + str(c))
        #print("Delta is flag: " + str(delta.is_flag_complex()))
        A = q_graphic_arrangement(delta,q)
        print(A)
        pd = int(compute_proj_dim(A,q))
        print("Its projective dimension is: "+str(pd))

        proj_dims.append(pd)
        tested_complexes.append(delta.facets())
        if (pd<=1 and w) or (not w and (pd>1 or (pd==1 and delta.dim()>1))):
            print("Strict PD bound holds")
        else: 
            if (pd==2 and w):
                print("only Loose PD bound holds")
                counter_examples.append((delta.facets(),"Only loose WC PD bound holds"))
            if (pd>2 and w):
                really_bad+=1
                print("neither PD bound holds!?")
                counter_examples.append((delta.facets(),"Neither WC PD bound holds!?"))
            if (pd==1 and not w and delta.dim()==1):
                really_bad+=1
                print("non weakly chordal 1-skeleton, but the strict PD bound holds!?")
                counter_examples.append((delta.facets(),"Non weakly chordal 1-skeleton, but the strict PD bound holds!?"))

        l = len(proj_dims)
        if (l % (t+2)==0) and include_clique and include_skeleta:
            skel = tested_complexes[l-t-2]
            clique = tested_complexes[l-1]
            pd_skel = proj_dims[l-t-2]
            pd_clique = proj_dims[l-1]
            if pd_skel+1==pd_clique:
                print("Clique bound is proper!")
                counter_examples.append((skel,clique,pd_skel,pd_clique,"Clique bound is proper!"))
            if pd_skel+1<pd_clique:
                print("Clique bound is too proper!?")
                counter_examples.append((skel,clique,pd_skel,pd_clique,"Clique bound is too proper!?"))
            if pd_skel<=pd_clique:
                print("Clique bound holds.")
            else:
                really_bad+=1
                print("Clique bound fails!?")
                counter_examples.append((skel,clique,pd_skel,pd_clique,"Clique bound fails!?"))

        progress += 1
        print("")

    #progress += 1
    print("One skeleton done!")
    print("")
#print(tested_complexes)
print(str(len(tested_complexes)) + " complexes on "+ str(n)+" vertices.")
print("Projective Dimensions: ")
print(proj_dims)
print("Counter Examples: ")
print(counter_examples)
print("Really Bad Counter Examples: " + str(really_bad))
