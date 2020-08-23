#include "cr_edgemidpoints.h"

template <typename DerivedV, typename DerivedF, typename DerivedE,
typename DerivedoE, typename Derivedmidpoints>
IGL_INLINE void
igl::cr_edgemidpoints(
                   const Eigen::MatrixBase<DerivedV>& V,
                   const Eigen::MatrixBase<DerivedF>& F,
                   const Eigen::MatrixBase<DerivedE>& E,
                   const Eigen::MatrixBase<DerivedoE>& oE,
                   Eigen::PlainObjectBase<Derivedmidpoints>& midpoints)
{
    assert(V.rows()==F.maxCoeff()+1 && "V does not seem to correspond to F.");
    
    const int m = E.maxCoeff()+1;
    
    midpoints.resize(m, 3);
    for(int i=0; i<F.rows(); ++i) {
        for(int j=0; j<3; ++j) {
            if(oE(i,j)<0)
                continue;
            midpoints.row(E(i,j)) = 0.5 * (V.row(F(i,(j+2)%3))
                                           + V.row(F(i,(j+1)%3)));
        }
    }
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::cr_edgemidpoints<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
