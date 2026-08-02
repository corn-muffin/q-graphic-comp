#This will print all RREFs of the orthogonal complement of itersections of hyperplanes
#This code scales very poorly

import itertools

def q_graphic_arrangement_norm(delta,q):
    n = len(delta.vertices())
    K = GF(q)
    hyperplanes = []

    for face in delta.face_iterator():
        if dim(face) == -1: #don't need to do anything for the empty set
            continue

        m = dim(face) + 1 # Number of vertices of face
        nonzero_K = [a for a in K if a != 0]
        V = [vector(K, l)  for l in itertools.product(nonzero_K, repeat=(m-1))] #all (m-1)-long lists of nonzero coefficients
        
        for v in V:
            v = vector([1] + list(v))
            j = 0
            a = [] # Coefficients

            # Insert coefficients from vector to entries of a indexed by face
            for i in [0..n-1]:
                if i in face:
                    a.append(v[j])
                    j += 1
                else:
                    a.append(0)
            hyperplanes.append(vector(a))
    return hyperplanes

def index_of_RREF(R):
     return (len(R.rows()),) + tuple(entry for row in R.rows() for entry in row)

def all_RREFs(delta,q):
    rrefs = set([])
    
    hyperplanes = q_graphic_arrangement_norm(delta,q)
    for A in Subsets([tuple(v) for v in hyperplanes]):
        M = matrix(vector(v) for v in A)
        M = M.echelon_form()
        M = M[:M.rank(), :]
        M.set_immutable()
        rrefs.add(M)
    return(rrefs)



#######example#########
q=3
rrefs = all_RREFs(SimplicialComplex([[0,1,2]]),q)
for R in sorted(rrefs, key=lambda M: index_of_RREF(M)):
    print(R,"\n")
