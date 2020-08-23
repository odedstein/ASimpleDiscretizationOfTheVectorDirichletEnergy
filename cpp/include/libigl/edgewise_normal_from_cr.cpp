#include "edgewise_normal_from_cr.h"

#include <igl/per_face_normals.h>
#include <igl/orient_halfedges.h>


template <typename DerivedN, typename DerivedE, typename DerivededgeN>
IGL_INLINE void
igl::edgewise_normal_from_cr(
                          const Eigen::MatrixBase<DerivedN>& N,
                          const Eigen::MatrixBase<DerivedE>& E,
                          Eigen::PlainObjectBase<DerivededgeN>& edgeN)
{
    assert(N.rows() == E.rows());
    
    edgeN.resize(E.maxCoeff()+1, 3);
    edgeN.setZero();
    for(int i=0; i<E.rows(); ++i)
        for(int j=0; j<3; ++j)
            edgeN.row(E(i,j)) += N.row(i);
    edgeN.rowwise().normalize();
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::edgewise_normal_from_cr<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
