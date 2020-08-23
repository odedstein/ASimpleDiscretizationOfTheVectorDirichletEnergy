// Std lib
#include <iostream>

// Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Project
#include <igl/cr_shape_operator.h>

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
    igl::matlab::mexErrMsgTxt(nlhs<=3,
                              "At most three output arguments are allowed.");
    
    igl::matlab::mexErrMsgTxt(nrhs >= 2,
                              "The number of input arguments must be at least 2.");
    igl::matlab::mexErrMsgTxt(nrhs <= 4,
                              "The number of input arguments must be at most 4.");
    igl::matlab::mexErrMsgTxt((nrhs == 2) || (nrhs == 4),
                              "Either E and oE are both provided, or none");
    const bool edgeLabelingProvided = (nrhs == 4);

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
    if(edgeLabelingProvided) {
        igl::matlab::mexErrMsgTxt(mxGetN(prhs[2]) == 3,
                                  "E has 3 columns");
        igl::matlab::mexErrMsgTxt(mxGetN(prhs[3]) == 3,
                                  "oE has 3 columns");
        igl::matlab::mexErrMsgTxt(mxGetM(prhs[2]) == mxGetM(prhs[1]),
                                  "E must have as many rows as F");
        igl::matlab::mexErrMsgTxt(mxGetM(prhs[3]) == mxGetM(prhs[1]),
                                  "oE must have as many rows as F");

        igl::matlab::parse_rhs_index(prhs+2,E);
        igl::matlab::parse_rhs_double(prhs+3,oE);
    }
    
    Eigen::MatrixXd S;
    
    //Construct matrices
    igl::cr_shape_operator(V, F, E, oE, S);
    
    //Output matrices
    igl::matlab::prepare_lhs_double(S, plhs+0);
    if(nlhs>=2)
        igl::matlab::prepare_lhs_index(E, plhs+1);
    if(nlhs>=3)
        igl::matlab::prepare_lhs_double(oE, plhs+2);
}
