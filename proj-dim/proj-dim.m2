-- Tools to compute the projective dimension of
-- different simplices with the same 1-skeleton.

needsPackage "Graphs"
needsPackage "SimplicialComplexes"
needsPackage "HyperplaneArrangements"
needsPackage "EdgeIdeals"

deleteFace = (C, ringElement) -> (

)

-- n: number of vertices
-- m: number of edges
-- l: number of random simplicial 
--   complexes you'd like too generate.
randomSimplices = (n, m, l) -> (
  if n<1 or m<0 or m> binomial(n,2) then (
    print "Your numerics are off!"
    return -1)
  R = QQ[apply(1..n, i->x_i)]
  G = randomGraph(R,m)
  C = cliqueComplex(G)
  d = dim(G)
  F = faces(C)
