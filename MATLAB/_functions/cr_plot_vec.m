function [vec,edgemps] = cr_plot_vec(V,F,E,oE,z)
%CR_PLOT_VEC Given a CR mesh structure V, F, E, oE, and a function z in the
%CR basis, return the embedded vectors vec at edgemidpoints edgemps

m = max(E(:));
aE = all_edges(F);
posEdges = find(oE(:)>0);
eE = nan(m,2); eE(E(posEdges),:) = aE(posEdges,:);

edgemps = 0.5 * (V(eE(:,1),:) + V(eE(:,2),:));

% Construct vecs to plot
edgeN = indexed_per_edge_normals(V, F, E);
vec = zeros(m,3);
for i=1:m
    edge = V(eE(i,2),:) - V(eE(i,1),:);
    edge = edge / norm(edge);
    edgePerp = edge * axisangle_to_rotmatrix(edgeN(i,:),pi/2)';
    vec(i,:) = edge*z(i) + edgePerp*z(i+m);
end

end

