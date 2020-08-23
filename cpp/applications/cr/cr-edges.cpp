// Std lib
#include <iostream>

// Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Project
#include <igl/doublearea.h>
#include <igl/cr_edges.h>
#include <igl/edgewise_normal_from_cr.h>
#include <igl/cr_edgemidpoints.h>

// MATLAB
#include <mex.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/MexStream.h>
#include <igl/C_STR.h>


void mexFunction(
                 int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    igl::matlab::mexErrMsgTxt(nlhs>=1,
                              "At least one output argument is needed.");
    igl::matlab::mexErrMsgTxt(nlhs<=2,
                              "At most two output arguments are possible.");
    
    igl::matlab::mexErrMsgTxt(nrhs == 4,
                              "4 input arguments are needed.");

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    const int dim = mxGetN(prhs[0]);
    igl::matlab::mexErrMsgTxt(dim == 3 || dim == 2,
                              "Mesh vertex list must be #V by 2 or 3 list of vertex positions");
    igl::matlab::mexErrMsgTxt(mxGetN(prhs[1]) == 3,
                              "Faces must have three vertices");
    igl::matlab::parse_rhs_double(prhs+0,V);
    igl::matlab::parse_rhs_index(prhs+1,F);
    
    Eigen::MatrixXi E, oE;
    igl::matlab::mexErrMsgTxt(mxGetN(prhs[2]) == 1,
                              "E is a column vector");
    igl::matlab::mexErrMsgTxt(mxGetN(prhs[3]) == 1,
                              "oE is a column vector");
    igl::matlab::mexErrMsgTxt(mxGetM(prhs[2]) == 3*mxGetM(prhs[1]),
                              "E must have as many rows as F halfedges");
    igl::matlab::mexErrMsgTxt(mxGetM(prhs[3]) == 3*mxGetM(prhs[1]),
                              "oE must have as many rows as F halfedges");
    igl::matlab::parse_rhs_index(prhs+2,E);
    igl::matlab::parse_rhs_double(prhs+3,oE);
    
    Eigen::MatrixXd faceN, edgeN, edges;
    
    //Construct edges
    igl::doublearea(V, F, faceN);
    igl::edgewise_normal_from_cr(faceN, E, edgeN);
    igl::cr_edges(V, F, E, oE, edgeN, edges);
    
    //Output matrices
    igl::matlab::prepare_lhs_double(edges, plhs);

     if(nlhs>=2) {
        //Construct and output midpoints
        Eigen::MatrixXd midpts;
        igl::cr_edgemidpoints(V, F, E, oE, midpts);
        igl::matlab::prepare_lhs_double(midpts, plhs+1);
     }
}