function E = all_edges(F)
%ALL_EDGES For a mesh F, return the from-to vertex index of all undirected
% edges with the following convention: first the edges [2,3], then the
% edges [3,1], then the edges [1,2]

E = [F(:,2) F(:,3); F(:,3) F(:,1); F(:,1) F(:,2)];

end

