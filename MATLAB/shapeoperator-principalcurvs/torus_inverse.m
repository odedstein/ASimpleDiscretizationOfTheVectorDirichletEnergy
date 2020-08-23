function [u,v] = torus_inverse(R, r, V)
%TORUS_INVERSE Given points V on a torus with large radius R and small
% radius r, return parametrization coordinates u, v

% Polar coordinates rho and param coord u
[u, rho] = cart2pol(V(:,1), V(:,2));

% atan2 to get coord v
v = atan2(V(:,3), rho-R);


end

