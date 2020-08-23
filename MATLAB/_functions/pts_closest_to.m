function inds = pts_closest_to(pts, X)
%PTS_CLOSEST_TO Finds the point in pts closest to each row of X and returns
%them as an index into pts

inds = nan(size(X,1),1);
for i=1:size(X,1)
    dists = normrow(pts - X(i,:));
    [~,minInd] = min(dists);
    inds(i) = minInd;
end

end

