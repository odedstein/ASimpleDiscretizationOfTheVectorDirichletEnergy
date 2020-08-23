function bd = bdry_edge_indices(E)
%BDRY_EDGE_INDICES Given halfedge indexing
%E, find the indices of all unique edges at the boundary.

%boundary edges occur only once in E, so find all these
bd = [];
for i=1:max(E(:))
    if(sum(E(:)==i)==1)
        bd = [bd;i];
    end
end

end

