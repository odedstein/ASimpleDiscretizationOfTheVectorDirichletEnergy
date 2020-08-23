% function [H, A, E, oE] = cr_holomorphic_energies_laplacian(V,F,E,oE)
% CR_VECTOR_LAPLACIAN Compute the Crouzeix-Raviart holomorphic and
% antiholomorphic energy matrix for a triangle mesh with vertices V and faces F.
% The arguments E and oE are optional. Both must be provided or none.
%  If they are provided, the method will not compute its own edge indexing,
%  and will use the edge indexing provided by E, oE.
%
% H 2*#edges by 2*#edges sparse matrix
%  the Crouzeix-Raviart holomorphic energy, two DOF per edge (first come all
%  the parallel ones, then all the normal)
% A 2*#edges by 2*#edges sparse matrix
%  the Crouzeix-Raviart antiholomorphic energy, two DOF per edge (first come all
%  the parallel ones, then all the normal)
% E #halfedges by 1 matrix
%  an edge convention that maps each unique halfedge to an edge index used
%  for indexing in L.
%  The halfedge indexing convention successively indexes all first halfedges,
%  then all second halfedges, then all third halfedges.
% oE #halfedges by 1 matrix
%  indicates the orientation of each edge
%
