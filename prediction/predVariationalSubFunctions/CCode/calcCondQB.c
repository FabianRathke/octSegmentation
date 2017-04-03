#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
int min(int A, int B) {
    if (A < B) {
        return A;
    } else {
        return B;
    }
}

int max(int A, int B) {
    if (A > B) {
        return A;
    } else {
        return B;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	/* Input variables */
	double *mu = mxGetPr(prhs[0]);
	double *var = mxGetPr(prhs[1]);
	int numRows = (int) mxGetScalar(prhs[2]);
	int numBounds = (int) mxGetScalar(prhs[3]);
	int numColumnsShape = (int) mxGetScalar(prhs[4]);
	int numColumnsPred = (int) mxGetScalar(prhs[5]);
    int *colA = (int*) mxGetData(prhs[6]);
    double *colAFac = mxGetPr(prhs[7]);
    int *colB = (int*) mxGetData(prhs[8]);
    double *colBFac = mxGetPr(prhs[9]);
  	double eps = mxGetScalar(prhs[10]);

	/* intern variables and pointers */
	float* condQB = NULL;
	int* cmin = NULL;
    int* cmax = NULL;	
	float factor,varinv,varAvg,muAvg;
	int i,j,k,idx,idx2,numNotZero,muFloor,startVal,stopVal,idxA,idxB;

  	/* 2-D matrix with [numBounds,numColumnsPred] */
	plhs[0] = mxCreateNumericMatrix(numRows,numBounds*numColumnsPred,mxSINGLE_CLASS,mxREAL);
    condQB = (float *) mxGetPr(plhs[0]);
	plhs[1] = mxCreateNumericMatrix(1,numBounds*numColumnsPred,mxINT32_CLASS,mxREAL);
    cmin =  (int *)mxGetData(plhs[1]);
	plhs[2] = mxCreateNumericMatrix(1,numBounds*numColumnsPred,mxINT32_CLASS,mxREAL);
    cmax =  (int *)mxGetData(plhs[2]);

	/* negative entropy of q_c */
	#pragma omp parallel for private(j,i,idx,idx2,startVal,stopVal,varinv,numNotZero,muFloor,factor,varAvg,muAvg,idxA,idxB)
	for (k=0; k < numBounds; k++) {
		for (j=0; j < numColumnsPred; j++) {
			/*idx = k*numColumnsShape + j;
			factor = 1/sqrt(2*3.1415926535897*var[idx]);
			muFloor = (int) mu[idx];*/
			idx = k*numColumnsPred + j;
			idxA = colA[j]+k*numColumnsShape; idxB = colB[j]+k*numColumnsShape;
			varAvg = (float) (colAFac[j]*var[idxA] + colBFac[j]*var[idxB]);
			muAvg = (float) (colAFac[j]*mu[idxA] + colBFac[j]*mu[idxB]);
			factor = 1/sqrtf(2*3.1415926535897*varAvg);
			muFloor = (int) muAvg;
			/* calculate rows for which the gaussian is larger than threshold */
			numNotZero = (int) ceil(abs(sqrt(-log(eps*factor)*2*varAvg)));
		  	startVal = max(muFloor-numNotZero,1);
			stopVal = min(muFloor+numNotZero,numRows);
			cmin[idx] = startVal-1;
			cmax[idx] = stopVal-1;
			idx2 = idx*numRows;
			varinv = -1/(2*varAvg);
			for (i=startVal; i <= stopVal; i++) {
				condQB[idx2 + i - 1] = factor*expf(varinv*(i-muAvg)*(i-muAvg));
			}
		}
	}
}
