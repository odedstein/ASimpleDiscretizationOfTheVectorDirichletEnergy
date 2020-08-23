% function [edges, midpts] = cr_edges(V,F,E,oE)
% CR_EDGES Compute a 3D representative for each edge, given a CR data structure.
% The arguments E and oE are optional. Both must be provided or none.
%  If they are provided, the method will not compute its own edge indexing,
%  and will use the edge indexing provided by E, oE.
%
% V, F: input mesh
% E #halfedges by 1 matrix
%  an edge convention that maps each unique halfedge to an edge index used
%  for indexing in L.
%  The halfedge indexing convention successively indexes all first halfedges,
%  then all second halfedges, then all third halfedges.
% oE #halfedges by 1 matrix
%  indicates the orientation of each edge
%
% edges: the embedded edge vectors (first the parallel, then the orthogonal
%  ones)
% midpoints: the 3D embedded midpoint of each edge
