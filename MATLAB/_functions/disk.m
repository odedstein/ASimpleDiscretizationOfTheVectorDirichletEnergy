function [V,F,b] = disk(n,ugly)
%DISK Constructs a mesh of a disk of radius 1 centered at 0 with n samples
%on the boundary
%If ugly==true, will generate a mesh with deformed triangles
%V,F: mesh
%b: boundary

if(nargin<2)
    ugly = false;
end

if(ugly)
    theta = linspace(0,2*pi,n+1)';
    theta = theta(1:end-1);
   
    a = n/6;
    b = 1;
    
    %rf = @(t) a*b ./ sqrt((b*cos(t)).^2 + (a*sin(t)).^2);
    %Vr = rf(theta).*[cos(theta) sin(theta)];
    Vr = [a*cos(theta) sin(theta)];
    Er = fliplr([1:size(Vr,1);2:size(Vr,1) 1]');
    
    flags = sprintf('-q30 -a%0.17f',(2*pi*sqrt(a*b)/n)^2);
    [V,F] = triangle(Vr,Er,[],'Flags',flags);
    
    %V = V ./ rf(atan2(V(:,2), V(:,1)));
    V(:,1) = V(:,1) / a;
    
    b = unique(outline(F));
    V(b,:) = V(b,:) ./ normrow(V(b,:));
else
    theta = linspace(0,2*pi,n+1)';
    theta = theta(1:end-1);

    Vr = [cos(theta) sin(theta)];
    Er = fliplr([1:size(Vr,1);2:size(Vr,1) 1]');
    
    flags = sprintf('-q30 -a%0.17f',(2*pi/n)^2);
    [V,F] = triangle(Vr,Er,[],'Flags',flags);
    b = unique(outline(F));
end

[V,I] = remove_unreferenced(V,F);
F = I(F);

end
