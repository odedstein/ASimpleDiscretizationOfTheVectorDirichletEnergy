
#include "cr_shape_operator.h"

#include "orient_halfedges.h"

#include "doublearea.h"
#include "squared_edge_lengths.h"
#include "per_face_normals.h"
#include "cr_embedding_derivative.h"

#include <Eigen/Geometry>

template <typename DerivedV, typename DerivedF, typename DerivedE,
typename DerivedOE, typename DerivedS>
IGL_INLINE void
igl::cr_shape_operator(
                  const Eigen::MatrixBase<DerivedV>& V,
                  const Eigen::MatrixBase<DerivedF>& F,
                  Eigen::PlainObjectBase<DerivedE>& E,
                  Eigen::PlainObjectBase<DerivedOE>& oE,
                  Eigen::PlainObjectBase<DerivedS>& S)
{
    if(E.rows()!=F.rows() || E.cols()!=F.cols() || oE.rows()!=F.rows() ||
       oE.rows()!=F.cols())
        orient_halfedges(F, E, oE);
    
    cr_shape_operator(V, F,
                      const_cast<const Eigen::PlainObjectBase<DerivedE>& >(E),
                      const_cast<const Eigen::PlainObjectBase<DerivedOE>& >(oE),
                      S);
}


template <typename DerivedV, typename DerivedF, typename DerivedE,
typename DerivedOE, typename DerivedS>
IGL_INLINE void
igl::cr_shape_operator(
                  const Eigen::MatrixBase<DerivedV>& V,
                  const Eigen::MatrixBase<DerivedF>& F,
                  const Eigen::MatrixBase<DerivedE>& E,
                  const Eigen::MatrixBase<DerivedOE>& oE,
                  Eigen::PlainObjectBase<DerivedS>& S)
{
    assert(V.cols()==3 && "Only works for surfaces embedded in 3D");
    assert(F.cols()==3 && "Faces have three vertices");
    assert(F.maxCoeff()+1 <= V.size() && "V too small for F");
    assert(E.rows()==F.rows() && E.cols()==F.cols() && oE.rows()==F.rows() &&
           oE.cols()==F.cols() && "Wrong dimension in edge vectors");
    
    using Mat = Eigen::Matrix<typename DerivedV::Scalar,
    Eigen::Dynamic, Eigen::Dynamic>;
    
    Mat l_sq, dA, N, dX;
    squared_edge_lengths(V, F, l_sq);
    doublearea(l_sq.array().sqrt().matrix(), dA);
    per_face_normals(V, F, N);
    cr_embedding_derivative(V, F, l_sq, dA, E, oE, dX);
    
    cr_shape_operator(V, F, dX, N, l_sq, dA, E, oE, S);
}


template <typename DeriveddX, typename DerivedN, typename DerivedV,
typename DerivedF, typename DerivedL_sq, typename DeriveddA,
typename DerivedE, typename DerivedOE, typename DerivedS>
IGL_INLINE void
igl::cr_shape_operator(
                  const Eigen::MatrixBase<DerivedV>& V,
                  const Eigen::MatrixBase<DerivedF>& F,
                  const Eigen::MatrixBase<DeriveddX>& dX,
                  const Eigen::MatrixBase<DerivedN>& N,
                  const Eigen::MatrixBase<DerivedL_sq>& l_sq,
                  const Eigen::MatrixBase<DeriveddA>& dA,
                  const Eigen::MatrixBase<DerivedE>& E,
                  const Eigen::MatrixBase<DerivedOE>& oE,
                  Eigen::PlainObjectBase<DerivedS>& S)
{
    assert(V.cols()==3 && "Only works for surfaces embedded in 3D");
    assert(F.cols()==3 && "Faces have three vertices");
    assert(F.maxCoeff()+1 <= V.size() && "V too small for F");
    assert(E.rows()==F.rows() && E.cols()==F.cols() && oE.rows()==F.rows() &&
           oE.cols()==F.cols() && "Wrong dimension in edge vectors");
    assert(l_sq.rows()==F.rows() && l_sq.cols()==3 && "l_sq dimensions wrong");
    assert(dA.size()==F.rows() && "dA dimensions wrong");
    assert(N.rows()==F.rows() && "N is a per-face normal vector");
    assert(N.cols()==3 && "Normal vectors must be embedded in 3D");
    
    using Scalar = typename DerivedS::Scalar;
    using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Vec3 = Eigen::Matrix<Scalar, 3, 1>;
    using Mat32 = Eigen::Matrix<Scalar, 3, 2>;
    using Mat3 = Eigen::Matrix<Scalar, 3, 3>;
    
    const int m = F.rows();
    const int nE = E.maxCoeff() + 1;
    
    assert(dX.rows()==2*nE && "dX has the number of CR basis fcts as cols.");
    assert(dX.cols()==3 && "dX has one derivative per coordinate fct");
    
    S.resize(m, 9);
    S.setZero();
    for(int f=0; f<m; ++f) {
        const Vec3& Nl = N.row(f);
        for(int e=0; e<3; ++e) {
            const Scalar eij=l_sq(f,e), ejk=l_sq(f,(e+1)%3),
            eki=l_sq(f,(e+2)%3); //These are squared quantities.
            const Scalar s_eij = sqrt(eij);
            const Scalar o = oE(f,e);
            const int i=F(f,(e+1)%3), j=F(f,(e+2)%3), k=F(f,e);
            
            // The vector eij and eij_perp as a 3x2 matrix of embedded 3d
            // vectors
            Mat32 vecs;
            vecs.col(0) = (V.row(j) - V.row(i)) / s_eij;
            vecs.col(1) = Eigen::AngleAxis<Scalar>(0.5*M_PI, Nl) * vecs.col(0);
            
            // Coordinate (1,0) of \nabla (dPhi_k)_{ij par} and (0,1) of
            // \nabla (dPhi_k)_{ij perp} in the vecs basis
            const Scalar x = o * 2. * s_eij / dA(f);
            
            // The simplified version of vecs * \nabla (dPhi_k)_{ij par} * vecs
            // and vecs * \nabla (dPhi_k)_{ij perp} * vecs, i.e. the nablas in
            // the embedded 3D basis
            const Mat3 Sl =
            // This first row is only ever nonzero on the boundary (and maybe
            // not even then?)
            vecs.col(1)*vecs.col(0).transpose() * Nl.dot(dX.row(E(f,e)))*x +
            vecs.col(1)*vecs.col(1).transpose() * Nl.dot(dX.row(E(f,e)+nE))*x;
            
            // Add to finished matrix
            S.row(f).template segment<3>(0) += Sl.row(0);
            S.row(f).template segment<3>(3) += Sl.row(1);
            S.row(f).template segment<3>(6) += Sl.row(2);
        }
    }
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation

#endif
