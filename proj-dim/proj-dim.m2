-- Tools to compute the projective dimension of
-- different simplices with the same 1-skeleton.

needsPackage "Graphs"
needsPackage "SimplicialComplexes"
needsPackage "HyperplaneArrangements"
needsPackage "EdgeIdeals"

deleteFace = (C, X) -> (
  --- if isMember(x,facets C) then (
  ---   return simplicalComplex(delete(x,facets C))
  --- );
  d= dim C;
  xInC := true;
  while xInC do (
    for Y in facets C do (
      if Y==X then (
        if #facets(C)=1 then (
          return simplicialComplex
        );
        return simplicalComplex(delete(Y,facets C))
      );
      if isMember(Y,ideal(X)) then (
        C = deleteFace(C,Y)
      );
    );
    xInC = false
  );
  return C;
)

randomComplex = (n,m,l

--- n: number of vertices
--- m: number of edges

randomGraph(n, m) -> (
  if n<1 or m<0 or m> binomial(n,2) then (
    print "Your numerics are off!"
    return -1)
  R = QQ[apply(1..n, i->x_i)]
  G = randomGraph(R,m)
)

randomSimplices

--- n: number of vertices
--- m: number of edges
--- p: probability of including any face
--- s: number of random simplicial 
---   complexes you'd like to
---   generate on the 1-skeleton.
randomSimplices = (n, m, p, s) -> (
  return randomSimplices(randomGraph(n,m),p,s)
)

n = 6 -- Number of vertices
d = 4 -- Dimension of complex
p = 0.5 -- Probability of including a d-face
s = 10 -- Number of simplicial complexes to generate

progress = 1
for delta in complexes:
    print("######################################### (" + str(progress) + "/" + str(s) + ")")


    progress += 1
    print("")
