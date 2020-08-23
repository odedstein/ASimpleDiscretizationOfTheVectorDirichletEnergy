function [para,perp] = project_onto_cr_space(vec,ind,V,F,E,oE)
%PROJECT_ONTO_CR_SPACE Project embedded vectors vec given at edge midpoints
%with indices ind into para and perp values corresponding to the cr
%structure given by V, F, E, oE

m = max(E(:));
aE = all_edges(F);
posEdges = find(oE(:)>0);
eE = nan(m,2); eE(E(posEdges),:) = aE(posEdges,:);

edgeN = indexed_per_edge_normals(V, F, E);
vedgetang = vec - dot(vec, edgeN(ind,:), 2).*edgeN(ind,:);
edge = V(eE(ind,2),:) - V(eE(ind,1),:);
edge = edge ./ normrow(edge);
edgePerp = nan(size(edge));
for i=1:size(edge,1)
    edgePerp(i,:) = edge(i,:) * ...
        (axisangle_to_rotmatrix(edgeN(ind(i),:),pi/2)');
end
para = dot(vedgetang,edge,2);
perp = dot(vedgetang,edgePerp,2);

end

