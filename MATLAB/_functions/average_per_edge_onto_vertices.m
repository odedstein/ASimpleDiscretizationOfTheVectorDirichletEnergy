function avgd = average_per_edge_onto_vertices(E,z)
%AVERADE_PER_EDGE_ONTO_VERTICES Given edges E, average the
%per-edge scalar data z (column vector) onto the vertices v

n = max(E(:));

edgesPerVec = nan(n,1);
for i=1:n
    edgesPerVec(i) = sum(E(:)==i);
end
avgd = full(sparse(E(:),ones(2*size(E,1),1), [z;z], n,1)) ./ edgesPerVec;

end

