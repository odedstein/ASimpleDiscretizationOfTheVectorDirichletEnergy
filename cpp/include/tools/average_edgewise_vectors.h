
#ifndef TOOLS_AVERAGE_EDGEWISE_VECTORS_H
#define TOOLS_AVERAGE_EDGEWISE_VECTORS_H

#include <Eigen/Core>

#include <igl/igl_inline.h>

namespace tools
    {
    
    // Averages edgewise vectors onto vertices
    //
    // Inputs:
    //  F      mesh on which to average
    //  E      map from each halfedge to edge index
    //  X      one vector per edge (same edge indexing as E)
    // Outputs:
    //  averaged    the averages vectors, per edge
    
    template <typename DerivedV, typename DerivedF, typename DerivedE,
    typename DerivedX, typename Derivedaveraged>
    IGL_INLINE void
    average_edgewise_vectors(
                         const Eigen::MatrixBase<DerivedV>& V,
                         const Eigen::MatrixBase<DerivedF>& F,
                         const Eigen::MatrixBase<DerivedE>& E,
                         const Eigen::MatrixBase<DerivedX>& X,
                         Eigen::PlainObjectBase<Derivedaveraged>& averaged);
    
    }

#ifndef IGL_STATIC_LIBRARY
#  include "average_edgewise_vectors.cpp"
#endif

#endif
