
#include "average_edgewise_vectors.h"

#include <igl/internal_angles.h>


template <typename DerivedV, typename DerivedF, typename DerivedE,
typename DerivedX, typename Derivedaveraged>
IGL_INLINE void
tools::average_edgewise_vectors(
                     const Eigen::MatrixBase<DerivedV>& V,
                     const Eigen::MatrixBase<DerivedF>& F,
                     const Eigen::MatrixBase<DerivedE>& E,
                     const Eigen::MatrixBase<DerivedX>& X,
                     Eigen::PlainObjectBase<Derivedaveraged>& averaged)
{
    using Scalar = typename DerivedV::Scalar;
    using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    
    const int n = V.rows();
    
    averaged.resize(n, 3);
    averaged.setZero();
    
    VectorX perVertAngs = VectorX::Zero(n);
    MatrixX angles;
    igl::internal_angles(V, F, angles);
    for(int i=0; i<F.rows(); ++i) {
        for(int j=0; j<3; ++j) {
            const int tail=F(i,(j+1)%3), head=F(i,(j+2)%3);
            const double angtail=angles(i,(j+1)%3),
            anghead=angles(i,(j+2)%3);
            perVertAngs(tail) += angtail;
            perVertAngs(head) += anghead;
            averaged.row(tail) += angtail*X.row(E(i,j));
            averaged.row(head) += anghead*X.row(E(i,j));
        }
    }
    averaged.array().colwise() /= perVertAngs.array();
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void tools::average_edgewise_vectors<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
