function N = torus_normal(R, r, V)
%TORUS_NORMAL Given points V on a torus with large radius R and small
% radius r, return the torus normal (embedded in 3D)

[u,v] = torus_inverse(R, r, V);
N = [cos(u) .* cos(v), ...
    sin(u) .* sin(v), ...
    sin(v)];

end

