function [V,F] = half_torus(n,m,s,l)
%TORUS Generate a half torus mesh with resolution parameters n and m controlling
%the density of the meshing, s,l controlling location of inner circle and inner radius;

theta = linspace(0,2*pi,n);
r = cos(linspace(pi,0,m))*l + s;

[R,T] = meshgrid(r,theta);
%%Create top half
Z_top = sqrt(1 - ((R-s)/l).^2)*l;
%%Convert to surface
[X,Y,Z] = pol2cart(T,R,Z_top);
%%Convert to mesh
[V,F] = surf_to_mesh(X,Y,Z);

end

