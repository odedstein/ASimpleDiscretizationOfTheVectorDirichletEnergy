function R = axisangle_to_rotmatrix(u,x)
%AXIXANGLE_TO_ROTMATRIX Given a 3D unit vector u and an angle x in
%radians, compute the rotation matrix associated with axis u and angle x.

cost = cos(x(1));
sint = sin(x(1));
u = u / norm(u);
ux = u(1);
uy = u(2);
uz = u(3);

R = [cost+ux*ux*(1-cost), ...
    ux*uy*(1-cost)-uz*sint, ...
    ux*uz*(1-cost)+uy*sint; ...
    uy*ux*(1-cost)+uz*sint, ...
    cost+uy*uy*(1-cost), ...
    uy*uz*(1-cost)-ux*sint; ...
    uz*ux*(1-cost)-uy*sint, ...
    uz*uy*(1-cost)+ux*sint, ...
    cost+uz*uz*(1-cost)];
 
end

