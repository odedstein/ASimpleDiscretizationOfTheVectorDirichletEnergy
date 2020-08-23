#ifndef IGL_CR_EDGES_H
#define IGL_CR_EDGES_H

#include <Eigen/Core>

#include <igl/igl_inline.h>

namespace igl
{
    
    // Compute all 3D edges in a CR mesh structure with correct indexing
    //
    // Inputs:
    //  V, F  mesh
    //  E, oE  CROF edge identifier matrices
    //  edgeN  edgewise normals
    // Outputs:
    //  edges  matrix containing an explicit embedded 3D CR edges
    
    template <typename DerivedV, typename DerivedF, typename DerivedE,
    typename DerivedoE, typename DerivededgeN, typename Derivededges>
    IGL_INLINE void
    cr_edges(
               const Eigen::MatrixBase<DerivedV>& V,
               const Eigen::MatrixBase<DerivedF>& F,
               const Eigen::MatrixBase<DerivedE>& E,
               const Eigen::MatrixBase<DerivedoE>& oE,
               const Eigen::MatrixBase<DerivededgeN>& edgeN,
               Eigen::PlainObjectBase<Derivededges>& edges);
    
}

#ifndef IGL_STATIC_LIBRARY
#  include "cr_edges.cpp"
#endif

#endif
