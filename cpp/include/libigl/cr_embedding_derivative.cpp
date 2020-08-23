
#include "cr_embedding_derivative.h"

#include "orient_halfedges.h"


template <typename DerivedV, typename DerivedF, typename DerivedE,
typename DerivedOE, typename DeriveddX>
IGL_INLINE void
igl::cr_embedding_derivative(
                             const Eigen::MatrixBase<DerivedV>& V,
                             const Eigen::MatrixBase<DerivedF>& F,
                             Eigen::PlainObjectBase<DerivedE>& E,
                             Eigen::PlainObjectBase<DerivedOE>& oE,
                             Eigen::PlainObjectBase<DeriveddX>& dX)
{
    if(E.rows()!=F.rows() || E.cols()!=F.cols() || oE.rows()!=F.rows() ||
       oE.rows()!=F.cols())
        orient_halfedges(F, E, oE);
    
    cr_embedding_derivative(V, F,
                            const_cast<const
                            Eigen::PlainObjectBase<DerivedE>& >(E),
                            const_cast<const
                            Eigen::PlainObjectBase<DerivedOE>& >(oE),
                            dX);
}



template <typename DerivedV, typename DerivedF, typename DerivedE,
typename DerivedOE, typename DeriveddX>
IGL_INLINE void
igl::cr_embedding_derivative(
                        const Eigen::MatrixBase<DerivedV>& V,
                        const Eigen::MatrixBase<DerivedF>& F,
                        const Eigen::MatrixBase<DerivedE>& E,
                        const Eigen::MatrixBase<DerivedOE>& oE,
                        Eigen::PlainObjectBase<DeriveddX>& dX)
{
    assert(V.cols()==3 && "Only works for surfaces embedded in 3D");
    assert(F.cols()==3 && "Faces have three vertices");
    assert(F.maxCoeff()+1 <= V.size() && "V too small for F");
    assert(E.rows()==F.rows() && E.cols()==F.cols() && oE.rows()==F.rows() &&
           oE.cols()==F.cols() && "Wrong dimension in edge vectors");
    
    using Mat = Eigen::Matrix<typename DerivedV::Scalar,
    Eigen::Dynamic, Eigen::Dynamic>;
    
    Mat l_sq, dA;
    squared_edge_lengths(V, F, l_sq);
    doublearea(l_sq.array().sqrt().matrix(), dA);
           
    cr_embedding_derivative(V, F, l_sq, dA, E, oE, dX);
}


template <typename DerivedV, typename DerivedF, typename Derivedl_sq,
typename DerivedE, typename DeriveddA, typename DerivedOE,
typename DeriveddX>
IGL_INLINE void
igl::cr_embedding_derivative(
                        const Eigen::MatrixBase<DerivedV>& V,
                        const Eigen::MatrixBase<DerivedF>& F,
                        const Eigen::MatrixBase<Derivedl_sq>& l_sq,
                        const Eigen::MatrixBase<DeriveddA>& dA,
                        const Eigen::MatrixBase<DerivedE>& E,
                        const Eigen::MatrixBase<DerivedOE>& oE,
                        Eigen::PlainObjectBase<DeriveddX>& dX)
{
    assert(V.cols()==3 && "Only works for surfaces embedded in 3D");
    assert(F.cols()==3 && "Faces have three vertices");
    assert(F.maxCoeff()+1 <= V.size() && "V too small for F");
    assert(E.rows()==F.rows() && E.cols()==F.cols() && oE.rows()==F.rows() &&
           oE.cols()==F.cols() && "Wrong dimension in edge vectors");
    assert(l_sq.rows()==F.rows() && l_sq.cols()==3 && "l_sq dimensions wrong");
    assert(dA.size()==F.rows() && "dA dimensions wrong");
    
    using Scalar = typename DeriveddX::Scalar;
    using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    
    const int m = F.rows();
    const int nE = E.maxCoeff() + 1;
    
    dX.resize(2*nE, 3);
    dX.setZero();
    VecX M(nE);
    M.setZero();
    for(int f=0; f<m; ++f) {
        for(int e=0; e<3; ++e) {
            const Scalar eij=l_sq(f,e), ejk=l_sq(f,(e+1)%3),
            eki=l_sq(f,(e+2)%3); //These are squared quantities.
            const Scalar s_eij=sqrt(eij);
            const Scalar o = oE(f,e);
            const int i=F(f,(e+1)%3), j=F(f,(e+2)%3), k=F(f,e);
            
            // This next row is only ever nonzero on the bdry, and is ignored
            // by the shape operator.
            dX.row(E(f,e)) += o * 1./(6.*s_eij) * dA(f) * (V.row(j)-V.row(i));
            dX.row(E(f,e)+nE) += o * ( -(
            eij*(V.row(j)+V.row(i)) + //This is only ever nonzero on the bdry
             (eki-ejk)*(V.row(j)-V.row(i)))
            / (12.*s_eij)
            + s_eij/6. * V.row(k));
            
            /*if(o>0)
                dX.row(E(f,e)) += (V.row(j) - V.row(i)) / s_eij;
            dX.row(E(f,e)+nE) += 0.5*o*(V.row(k) - 0.5*(V.row(i)+V.row(j)))
            * s_eij / dA(f);*/
            
            M(E(f,e)) += dA(f) / 6.;
        }
    }
    
    for(int i=0; i<3; ++i) {
        dX.block(0,i,nE,1).array() /= M.array();
        dX.block(nE,i,nE,1).array() /= M.array();
    }
}



#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation

#endif
