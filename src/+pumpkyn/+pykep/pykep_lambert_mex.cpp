#include "mex.h"
#include "lambert_original.hpp"
#include "const.hpp"
#include "vector3d.hpp"
#include <stdio.h>
#include <string.h>
#include <iostream>

// To Comile Source Code in MATLAB:
// mex COPTIMFLAGS="-Ofast -fp:fast" pykep_lambert_mex.cpp elliptic_orbit.cpp lambert_original.cpp
//
// [v1,v2] = pykep_lambert_mex(r1,r2,tof,mu,dir,nmax);
//

/* the gateway function to lambert source code*/
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) 
{
    using namespace astro_cpp;
    std::ios::sync_with_stdio(false);
    std::cout.precision(16);
    double *r1, *r2, *v1, *v2,tof;
    astro_cpp::Vector3D r1In, r2In, v1Out, v2Out;
    double r1vec[3],r2vec[3],v1vec[3],v2vec[3];
    double mu,dir;
    int i, max_revs = 0;
    long int nmax;
    size_t mrows, ncols;
    bool is_retrograde = false;
    
    /*  create a pointer to required inputs */
    r1 = mxGetPr(prhs[0]);
    r2 = mxGetPr(prhs[1]);
   tof = mxGetScalar(prhs[2]);
    mu = mxGetScalar(prhs[3]);
   dir = mxGetScalar(prhs[4]);
  nmax = mxGetScalar(prhs[5]);
  
    /* Set the retrograde flag */
    if (dir < 0.0){is_retrograde = true;}
    
    // Input initial/finalpositions:
    r1In.x = r1[0];
    r1In.y = r1[1];
    r1In.z = r1[2];
    r2In.x = r2[0];
    r2In.y = r2[1];
    r2In.z = r2[2];
    
    // Create an instance of Lambert Problem
     lambert_original lambert(r1In,r2In, tof, mu, is_retrograde, nmax);
    
    // Pull out required array size:
    ncols =  lambert.get_Nmax()*2 + 1;
    mrows = 3;
    
     /*  set the output pointers to the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
    
    /*  create a C pointer to a copy of the output matrix */
        v1 = mxGetPr(plhs[0]);
        v2 = mxGetPr(plhs[1]); 
        
   // Pull out velocity solution(s) and assign to pre-allocated arrays:
   for (i=0; i<ncols; i++)
   {
            int cOffset = i*3;
                  v1Out = lambert.get_v1()[i]; 
                  v2Out = lambert.get_v2()[i];
              *(v1+cOffset+0) = v1Out.x;
              *(v1+cOffset+1) = v1Out.y;
              *(v1+cOffset+2) = v1Out.z;
              *(v2+cOffset+0) = v2Out.x;
              *(v2+cOffset+1) = v2Out.y;
              *(v2+cOffset+2) = v2Out.z;
   }     
}