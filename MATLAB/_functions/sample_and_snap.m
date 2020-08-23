function [inds,presentInds] = sample_and_snap(pts,V,F,res,present,whiteInstead)
%SAMPLE_AND_SNAP sample random points on the mesh V,F at resolution res,
%and then snap to points pts. Returns indices into pts.
%Optional: if present is set, will incorporate it into inds

A = 0.5 * sum(doublearea(V,F));
n = ceil(A*res);

if(nargin<6)
    whiteInstead = false;
end

coloru = 'blue';
if(whiteInstead)
    coloru = 'white';
end

N = random_points_on_mesh(V,F,n, 'Color',coloru, 'MaxIter',1000000);
if(nargin>4)
    presentInds = nan(size(present,1),1);
    for i=1:size(present,1)
        [~,minind] = min(normrow(N-present(i,:)));
        N(minind,:) = present(i,:);
        presentInds(i) = minind;
    end
end
inds = snap_points(N,pts);

end

