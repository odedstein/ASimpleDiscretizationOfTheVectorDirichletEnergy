function N = indexed_per_edge_normals(V,F,E)
%PER_EDGE_NORMALS Given a halfedges-to-edges map E, create per-edge normals
%with the same indexing

m = max(E(:));
dim = size(V,2);

faceN = normals(V,F);

N = zeros(m,dim);
for j=1:size(faceN,1)
    N(E(j,:),:) = N(E(j,:),:) + repmat(faceN(j,:), [3 1]);
end

N = N ./ normrow(N);

end

