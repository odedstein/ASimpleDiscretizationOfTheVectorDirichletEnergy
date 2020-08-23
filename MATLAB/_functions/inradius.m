function r = inradius(V,F)
%INRADIUS Given a triangle mesh V,F, compute the inradius of each face

aE = V(F(:,3),:) - V(F(:,2),:);
a = normrow(aE);
bE = V(F(:,1),:) - V(F(:,3),:);
b = normrow(bE);
cE = V(F(:,2),:) - V(F(:,1),:);
c = normrow(cE);

r = 0.5 * sqrt(((b+c-a).*(c+a-b).*(a+b-c)) ./ (a+b+c));

end

