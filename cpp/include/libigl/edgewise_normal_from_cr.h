#ifndef TOOLS_EDGEWISE_NORMAL_FROM_CR_H
#define TOOLS_EDGEWISE_NORMAL_FROM_CR_H

#include <Eigen/Core>

#include <igl/igl_inline.h>

namespace igl
{
    
    // Compute the edge-wise normal on each edge
    //
    // Inputs:
    //  N  facewise normals
    //  E  CR edge identifier matrix
    // Outputs:
    //  edgeN  edgewise normals
    
    template <typename DerivedN, typename DerivedE, typename DerivededgeN>
    IGL_INLINE void
    edgewise_normal_from_cr(
                            const Eigen::MatrixBase<DerivedN>& N,
                            const Eigen::MatrixBase<DerivedE>& E,
                            Eigen::PlainObjectBase<DerivededgeN>& edgeN);
    
}

#ifndef IGL_STATIC_LIBRARY
#  include "edgewise_normal_from_cr.cpp"
#endif

#endif
