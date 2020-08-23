function avgd = vertexwise_plot_vec(N,X,z)
%vertexwise_plot_vec Given a per-vertex complex vector z, per-vertex normals N,
% and a basis for the tangent space X, return the embedded vectors vec at
% vertices

avgd = zeros(size(z,1),3);
for i=1:size(z,1)
    R = axisangle_to_rotmatrix(N(i,:), angle(z(i)));
    avgd(i,:) = abs(z(i))*(X(i,:)*(R'));
end

end

