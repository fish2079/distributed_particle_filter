/*
   Function call in matlab: replicatedIndices = systematicResamp(particles_weights,numParticles,uniformRndNum)
*/

#include "mex.h"

void mexFunction(int n_output, mxArray *output[], int n_input, const mxArray *input[])
{
    int m;
    int index = 0;
    int numParticles = (int) *mxGetPr(input[1]);
    
    int* replicatedIndices;

    double cumSum =  0.0; /*initializing cumulative sum of weights*/
    double oneOverNumPart = 1.0/numParticles;
    double uniformRndNum = (double) *mxGetPr(input[2]);

    double* particles_weights = mxGetPr(input[0]);
	
    /* Check the number of inputs and output arguments:*/
    if(n_output!=1) mexErrMsgTxt("Wrong number of output variables!");
	if(n_input!=3) mexErrMsgTxt("Wrong number of input variables!");
    
    uniformRndNum /= numParticles;
    
    /* Create output matrix:*/
    output[0] = mxCreateNumericArray(mxGetNumberOfDimensions(input[0]), mxGetDimensions(input[0]), mxINT32_CLASS, mxREAL);
    
    replicatedIndices = mxGetData(output[0]);
    
	/*------------------ Start of routine ---------------------------*/
    
    for(m=0;m<numParticles;m++)
    {   
        cumSum += particles_weights[m];
        while(cumSum > uniformRndNum)
        {
            uniformRndNum += oneOverNumPart;
            replicatedIndices[index] = m+1; /*m+1 -> due to Matlab indexing...*/
            index ++;
        }
    }
	    
	/* ------------------ End of routine -----------------------------*/
}
