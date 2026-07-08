## Code for computing the number of d dimensional subspaces intersecting W_j(q)'s trivially.

def avoids_arrangement(subspace,arrangement):
    for w in arrangement:
        if subspace.intersection(w).dimension() > 0:
            return False
    
    return True

def count_subspaces(d,n,q,facets):
    V = VectorSpace(GF(q),n)
    subspaces = V.subspaces(dim=d)

    W = []
    for sigma_j in facets:
        W.append(V.subspace([V.basis()[i] for i in sigma_j]))
    
    count = 0
    for H in subspaces:
        if avoids_arrangement(H,W):
            count += 1

    return count

Q = [2,3,5,7,11]

for q in Q:
    print(count_subspaces(2,5,q,[[0,1],[1,2,3],[3,4]]))