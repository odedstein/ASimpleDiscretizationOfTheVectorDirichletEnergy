function [V, F] = torus(R, r, n, m, varargin)
%TORUS Creates a torus mesh with large radius R, small radius r, n segments
% in the horizontal direction and m segments in the vertical direction.
% if the varargin argument is larger , will disturb coordinates a bit
% to give slightly irregular parametrization


%Setting up a torus mesh
u=linspace(0,2*pi,n);
v=linspace(0,2*pi,m);
[u,v]=meshgrid(u,v);

%perturb?
d = 0;
if(nargin>0)
    if(varargin{1}>0)
        d = varargin{1};
    end
end

if(d > 0)
    ru = ((2*pi)/n)^2 * 0.25 * d;
    rv = ((2*pi)/m)^2 * 0.25 * d;
    %ru = ((2*pi)/n) * 0.25 * d;
    %rv = ((2*pi)/m) * 0.25 * d;
    
    rng(34);
    randus = ru*rand(m,n);
    randus(2:2:end) = -randus(2:2:end);
    u = u+randus;
    u(:,end) = u(:,1);
    u(end,:) = u(1,:);
    randvs = rv*(rand(m,n)-0.5);
    v = v+randvs;
    v(:,end) = v(:,1);
    v(end,:) = v(1,:);
end


x=(R+r.*cos(v)).*cos(u);
y=(R+r.*cos(v)).*sin(u);
z=r.*sin(v);

[V, F] = surf_to_mesh(x, y, z);

end

