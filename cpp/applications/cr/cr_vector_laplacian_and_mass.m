% function [L, M, E, oE] = cr_vector_laplacian(V,F,E,oE)
% CR_VECTOR_LAPLACIAN Compute the Crouzeix-Raviart vector Laplacian and mass
%  matrix for a triangle mesh with vertices V and faces F.
% The arguments E and oE are optional. Both must be provided or none.
%  If they are provided, the method will not compute its own edge indexing,
%  and will use the edge indexing provided by E, oE.
%
% L 2*#edges by 2*#edges sparse matrix
%  the Crouzeix-Raviart vector Laplacian, two DOF per edge (first come all
%  the parallel ones, then all the normal)
% M 2*#edges by 2*#edges sparse matrix
%  the Crouzeix-Raviart vector mass matrix
% E #halfedges by 1 matrix
%  an edge convention that maps each unique halfedge to an edge index used
%  for indexing in L.
%  The halfedge indexing convention successively indexes all first halfedges,
%  then all second halfedges, then all third halfedges.
% oE #halfedges by 1 matrix
%  indicates the orientation of each edge
%
