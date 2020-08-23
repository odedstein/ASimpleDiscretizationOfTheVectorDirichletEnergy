#ifndef IGL_CR_EMBEDDING_DERIVATIVE_H
#define IGL_CR_EMBEDDING_DERIVATIVE_H

#include "igl_inline.h"

#include <Eigen/Core>
#include <Eigen/Sparse>


namespace igl
{
    // Computes the exterior derivative of the (coordinatewise) embedding
    //  function
    //
    // Inputs:
    //  V, F: input mesh
    //  E: a mapping from each halfedge to each edge.
    //  oE: the orientation of each halfedge compared to the orientation of the
    //      actual edge.
    //
    // Outputs:
    //  dX: the derivative of the mesh embedding function given as a three
    //      vector CR functions. This is a 2#nedges * 3 matrix
    //  E, oE: these are computed if they are not present.
    //

    template <typename DerivedV, typename DerivedF, typename DerivedE,
    typename DerivedOE, typename DeriveddX>
    IGL_INLINE void
    cr_embedding_derivative(
                            const Eigen::MatrixBase<DerivedV>& V,
                            const Eigen::MatrixBase<DerivedF>& F,
                            const Eigen::MatrixBase<DerivedE>& E,
                            const Eigen::MatrixBase<DerivedOE>& oE,
                            Eigen::PlainObjectBase<DeriveddX>& dX);

    template <typename DerivedV, typename DerivedF, typename DerivedE,
    typename DerivedOE, typename DeriveddX>
    IGL_INLINE void
    cr_embedding_derivative(
                            const Eigen::MatrixBase<DerivedV>& V,
                            const Eigen::MatrixBase<DerivedF>& F,
                            Eigen::PlainObjectBase<DerivedE>& E,
                            Eigen::PlainObjectBase<DerivedOE>& oE,
                            Eigen::PlainObjectBase<DeriveddX>& dX);


    // Version with precomputed inputs
    //
    // Inputs:
    //  V, F: input mesh
    //  l_sq: squared edge lengths of each halfedge
    //  dA: double area of each face
    //  E: a mapping from each halfedge to each edge.
    //  oE: the orientation of each halfedge compared to the orientation of the
    //      actual edge.
    //
    // Outputs:
    //  dX: the derivative of the mesh embedding function given as a three
    //      vector CR functions. This is a 2#nedges * 3 matrix
    //  E, oE: these are computed if they are not present.
    //

    template <typename DerivedV, typename DerivedF, typename Derivedl_sq,
    typename DerivedE, typename DeriveddA, typename DerivedOE,
    typename DeriveddX>
    IGL_INLINE void
    cr_embedding_derivative(
                            const Eigen::MatrixBase<DerivedV>& V,
                            const Eigen::MatrixBase<DerivedF>& F,
                            const Eigen::MatrixBase<Derivedl_sq>& l_sq,
                            const Eigen::MatrixBase<DeriveddA>& dA,
                            const Eigen::MatrixBase<DerivedE>& E,
                            const Eigen::MatrixBase<DerivedOE>& oE,
                            Eigen::PlainObjectBase<DeriveddX>& dX);
    
    
    
}


#ifndef IGL_STATIC_LIBRARY
#  include "cr_embedding_derivative.cpp"
#endif

#endif

