#ifndef TOOLS_CR_EDGEMIDPOINTS_H
#define TOOLS_CR_EDGEMIDPOINTS_H

#include <Eigen/Core>

#include <igl/igl_inline.h>

namespace igl
    {
    
    // Compute all edges in a CROF mesh structure with correct indexing
    //
    // Inputs:
    //  V, F  mesh
    //  E, oE  CROF edge identifier matrices
    // Outputs:
    //  midpoints  matrix containing the explicit embedded 3D CROF edges
    
    template <typename DerivedV, typename DerivedF, typename DerivedE,
    typename DerivedoE, typename Derivedmidpoints>
    IGL_INLINE void
    cr_edgemidpoints(
                       const Eigen::MatrixBase<DerivedV>& V,
                       const Eigen::MatrixBase<DerivedF>& F,
                       const Eigen::MatrixBase<DerivedE>& E,
                       const Eigen::MatrixBase<DerivedoE>& oE,
                       Eigen::PlainObjectBase<Derivedmidpoints>& midpoints);
    
    }

#ifndef IGL_STATIC_LIBRARY
#  include "cr_edgemidpoints.cpp"
#endif

#endif
