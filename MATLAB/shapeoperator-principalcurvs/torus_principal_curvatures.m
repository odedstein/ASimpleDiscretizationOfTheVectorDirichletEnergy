function [k1, k2, u1, u2] = torus_principal_curvatures(R, r, V)
%TORUS_PRINCIPAL_CURVATURES Given points V on a torus with large radius R
% and small radius r, return the torus principal curvatures k1, k2
% and the corresponding curvature directions u1, u2

[u,v] = torus_inverse(R, r, V);
n = size(u,1);

ks = [ cos(v) ./ (R + r*cos(v)); ...
    1/r * ones(n,1)];
us = [-sin(u).*cos(v), cos(u).*cos(v), zeros(n,1); ...
    -cos(u).*sin(v), -sin(u).*sin(v), cos(v)];
us = us ./ normrow(us);

%Sort by smallest / largest curvature
I1 = (1:n)';
I2 = ((n+1):(2*n))';
uorv = abs(ks(I1)) > abs(ks(I2));

k1 = uorv.*ks(I2) + (~uorv).*ks(I1);
k2 = (~uorv).*ks(I2) + uorv.*ks(I1);
u1 = uorv.*us(I2,:) + (~uorv).*us(I1,:);
u2 = (~uorv).*us(I2,:) + uorv.*us(I1,:);

end

