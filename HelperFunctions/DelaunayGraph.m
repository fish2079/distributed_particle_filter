function A = DelaunayGraph(X)

DT = delaunayTriangulation(X);
E = DT.edges;
A = sparse([E(:,1); E(:,2)], [E(:,2); E(:,1)], ones(2*size(E,1), 1));
A = full(A);