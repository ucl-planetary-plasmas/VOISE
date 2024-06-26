/**
 * @file
 * @brief This is a MEX function. It should only be compiled by the compileMEX.m
 * Matlab script. Removes multiple seeds from Voronoi diagram.
 *
 * @date Created 10/07/18
 * @author Jack Scantlebury
 */

#include "mexIncludes.h"

/**
 * @defgroup removeSeedFromVDBatch removeSeedFromVDBatch
 * @brief Removes multiple seeds from Voronoi diagram.
 *
 * This function interfaces with Matlab, taking the VD struct and performing
 * the removeSeed function repeatedly to remove seeds from the VD before giving
 * the result back to Matlab. The seeds to be removed are supplied as a list of
 * seed IDs.
 *
 * This is a MEX function. As such, the inputs and outputs are constricted to
 * the following:
 *
 * - nlhs: Number of outputs
 *
 * - plhs: Pointer to outputs
 *
 * - nrhs: Number of inputs
 *
 * - prhs: Pointer to inputs
 *
 * In Matlab, this corresponds to the following parameters and outputs:
 * @param VD Voronoi diagram struct
 * @param S List of seed IDs to be removed
 * @returns void
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

  // Sanity checks on inputs
  if (nlhs != 1 || nrhs != 2) {
    mexErrMsgTxt(
      " Invalid number of input and output arguments");
    return;
  }
  real *SkArray = mxGetDoubles(prhs[1]);
  int nSeeds = mxGetNumberOfElements(prhs[1]);

  // Grab VD data from ML struct
  vd outputVD = grabVD(prhs, 0);

  // Add seeds to VD
  for(auto i=0; i<nSeeds; ++i) {
    try {
      removeSeed(outputVD, SkArray[i]);
    }
    catch (SKIZException &e) {
      mexErrMsgTxt(e.what());
    }
  }

  // Push modified VD to ML VD struct (handles memory allocation etc)
  pushVD(outputVD, plhs);
}
