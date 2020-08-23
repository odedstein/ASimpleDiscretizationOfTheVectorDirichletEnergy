#ifndef IGL_CR_SHAPE_OPERATOR_H
#define IGL_CR_SHAPE_OPERATOR_H

#include "igl_inline.h"

#include <Eigen/Core>
#include <Eigen/Sparse>


namespace igl
{
    // Computes the CR shape operator.
    // This shape operator is constant per-face, and is returned as a 3x3 matrix
    //  representing the not-integrated shape operator.
    //
    // Inputs:
    //  V, F: input mesh
    //  E: a mapping from each halfedge to each edge.
    //  oE: the orientation of each halfedge compared to the orientation of the
    //      actual edge.
    //
    // Outputs:
    //  S: a 3x3-matrix Sl for each face, such that X'*Sl*Y = ShapeOp(X,Y) for
    //     each (embedded in 3D) tangent vector X,Y, , and N'*Sl*N = 0 for the
    //     normal vector N.
    //     S is a #nfaces * 9 matrix, where S(f,:) = [Sl(1,:) Sl(2,:) Sl(3,:)]
    //     for each face f and per-face matrix Sl.
    //  E, oE: these are computed if they are not present.

    template <typename DerivedV, typename DerivedF, typename DerivedE,
    typename DerivedOE, typename DerivedS>
    IGL_INLINE void
    cr_shape_operator(
                      const Eigen::MatrixBase<DerivedV>& V,
                      const Eigen::MatrixBase<DerivedF>& F,
                      const Eigen::MatrixBase<DerivedE>& E,
                      const Eigen::MatrixBase<DerivedOE>& oE,
                      Eigen::PlainObjectBase<DerivedS>& S);

    template <typename DerivedV, typename DerivedF, typename DerivedE,
    typename DerivedOE, typename DerivedS>
    IGL_INLINE void
    cr_shape_operator(
                      const Eigen::MatrixBase<DerivedV>& V,
                      const Eigen::MatrixBase<DerivedF>& F,
                      Eigen::PlainObjectBase<DerivedE>& E,
                      Eigen::PlainObjectBase<DerivedOE>& oE,
                      Eigen::PlainObjectBase<DerivedS>& S);




    // Version that uses some precomputed inputs
    //  dX: the derivative of the mesh embedding function given as a three
    //      vector CR functions. This is a 2#nedges * 3 matrix
    //  N: the per-face normal
    //  F: input mesh connectivity
    //  l_sq: squared edge lengths of each halfedge
    //  dA: double area of each face
    //  E: a mapping from each halfedge to each edge.
    //  oE: the orientation of each halfedge compared to the orientation of the
    //      actual edge.
    //
    // Outputs:
    //  S: a 3x3-matrix Sl for each face, such that X'*Sl*Y = ShapeOp(X,Y) for
    //     each (embedded in 3D) tangent vector X,Y, , and N'*Sl*N = 0 for the
    //     normal vector N.
    //     S is a #nfaces * 9 matrix, where S(f,:) = [Sl(1,:) Sl(2,:) Sl(3,:)]
    //     for each face f and per-face matrix Sl.

    template <typename DeriveddX, typename DerivedN, typename DerivedV,
    typename DerivedF, typename DerivedL_sq, typename DeriveddA,
    typename DerivedE, typename DerivedOE, typename DerivedS>
    IGL_INLINE void
    cr_shape_operator(
                      const Eigen::MatrixBase<DerivedV>& V,
                      const Eigen::MatrixBase<DerivedF>& F,
                      const Eigen::MatrixBase<DeriveddX>& dX,
                      const Eigen::MatrixBase<DerivedN>& N,
                      const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                      const Eigen::MatrixBase<DeriveddA>& dA,
                      const Eigen::MatrixBase<DerivedE>& E,
                      const Eigen::MatrixBase<DerivedOE>& oE,
                      Eigen::PlainObjectBase<DerivedS>& S);

}


#ifndef IGL_STATIC_LIBRARY
#  include "cr_shape_operator.cpp"
#endif

#endif

