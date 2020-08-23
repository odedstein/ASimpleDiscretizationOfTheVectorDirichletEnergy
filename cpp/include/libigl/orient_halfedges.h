#ifndef IGL_ORIENT_HALFEDGES_H
#define IGL_ORIENT_HALFEDGES_H

#include "igl_inline.h"
#include <Eigen/Core>


namespace igl
{
    // Orients halfedges for a mesh
    //
    // Inputs:
    //  F: input mesh connectivity
    //
    // Outputs:
    //  E: a mapping from each halfedge to each edge.
    //  oE: the orientation of each halfedge compared to the orientation of the
    //      actual edge.
    
    template <typename DerivedF, typename DerivedE, typename DerivedOE>
    IGL_INLINE void
    orient_halfedges(
                    const Eigen::MatrixBase<DerivedF>& F,
                    Eigen::PlainObjectBase<DerivedE>& E,
                    Eigen::PlainObjectBase<DerivedOE>& oE);
    
}
    
    
#ifndef IGL_STATIC_LIBRARY
#  include "orient_halfedges.cpp"
#endif

#endif
