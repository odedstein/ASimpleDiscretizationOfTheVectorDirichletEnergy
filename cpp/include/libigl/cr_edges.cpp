#include "cr_edges.h"

#include <Eigen/Geometry>

#include <iostream>

template <typename DerivedV, typename DerivedF, typename DerivedE,
typename DerivedoE, typename DerivededgeN, typename Derivededges>
IGL_INLINE void
igl::cr_edges(
           const Eigen::MatrixBase<DerivedV>& V,
           const Eigen::MatrixBase<DerivedF>& F,
           const Eigen::MatrixBase<DerivedE>& E,
           const Eigen::MatrixBase<DerivedoE>& oE,
           const Eigen::MatrixBase<DerivededgeN>& edgeN,
           Eigen::PlainObjectBase<Derivededges>& edges)
{
    edges.resize(2*edgeN.rows(), 3);
    edges.setZero();
    
    for(int i=0; i<F.rows(); ++i) {
        for(int j=0; j<3; ++j) {
            if(oE(i,j)<0)
                continue;
            const Eigen::AngleAxisd rotmat(0.5*M_PI, edgeN.row(E(i,j)));
            edges.row(E(i,j)) = (V.row(F(i,(j+2)%3))
                                 - V.row(F(i,(j+1)%3))).normalized();
            edges.row(E(i,j)+edgeN.rows()) = rotmat*edges.row(E(i,j)).transpose();
        }
    }
    
#ifndef NDEBUG
    for(int i=0; i<edges.rows(); ++i)
        assert(edges.row(i).norm() > 0.5);
#endif
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::cr_edges<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
