#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	/* Input variables */
	double *singleton = mxGetPr(prhs[0]);
	double *prediction = mxGetPr(prhs[1]);
	
	/* intern variables and pointers */
	double* q_c_singleton = NULL;
	double* q_c_data = NULL;
	double tmp;

	const mwSize *dim_array;
	dim_array = mxGetDimensions(prhs[0]);

	int i,j,k,idx;
	int numRows = dim_array[0];
	int numBounds = dim_array[1];
	int numColumnsPred = dim_array[2];

  	/* 2-D matrix with [numBounds,numColumnsPred] */
	plhs[0] = mxCreateDoubleMatrix(1,numBounds*numColumnsPred,mxREAL);
    q_c_singleton = mxGetPr(plhs[0]);
	/* 2-D matrix with [numColumnsPred,numBounds] */
  	plhs[1] = mxCreateDoubleMatrix(1,numBounds*numColumnsPred,mxREAL);
    q_c_data = mxGetPr(plhs[1]);


	/* negative entropy of q_c */
	#pragma omp parallel for private(tmp,k,i,idx)
	for (j=0; j < numColumnsPred; j++) {
        for (k=0; k < numBounds; k++) { 
			tmp = 0;
			idx = j*numBounds*numRows + k*numRows;
			for (i=0; i < numRows; i++) {
				if (singleton[idx + i] > 0) {
					tmp += singleton[idx+i]*log(singleton[idx+i]);
				}
			}
			q_c_singleton[k + j*numBounds] = -tmp;
		} 
	}

	/* data term */
	#pragma omp parallel for private(tmp,k,i,idx)
	for (j=0; j < numColumnsPred; j++) {
/*		idx = j*numBounds*numRows;
		tmp = 0;
		for (i=0; i < numRows; i++) {
			if (singleton[idx + i] > 0 & prediction[idx+i] > 0) {
				tmp += singleton[idx + i]*log(prediction[idx + i]);
			}
		}
        q_c_data[j] = -tmp;*/
		for (k=0; k < numBounds; k++) {
			for (j=0; j < numColumnsPred; j++) {
				idx = j*numBounds*numRows + k*numRows;
				tmp = 0;
				for (i=0; i < numRows; i++) {
					if (singleton[idx + i] > 0 && prediction[idx+i] > 0) {
						tmp += singleton[idx + i]*log(prediction[idx + i]);
					}
				}
				q_c_data[k*numColumnsPred+j] = -tmp;
			}
		}
	}
}
