% function [S, E, oE] = cr_shape_operator(V,F,E,oE)
% CR_SHAPE_OPERATOR Compute the pointwise per-face (constant) CR shape operator
%  of the mesh V,F.
% The arguments E and oE are optional. Both must be provided or none.
%  If they are provided, the method will not compute its own edge indexing,
%  and will use the edge indexing provided by E, oE.

% S #faces by 9 matrix
%  the shape operator.
%  Each row of S is a 3x3 matrix corresponding to the local shape operator Sl
%  for that face, S(f,:) = [Sl(1,:), Sl(2,:), Sl(3,:)]
% E #halfedges by 1 matrix
%  an edge convention that maps each unique halfedge to an edge index used
%  for indexing in L.
%  The halfedge indexing convention successively indexes all first halfedges,
%  then all second halfedges, then all third halfedges.
% oE #halfedges by 1 matrix
%  indicates the orientation of each edge
%
