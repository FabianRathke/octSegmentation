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
	double *prec = mxGetPr(prhs[1]);
	int numRows = (int) mxGetScalar(prhs[2]);
	int numBounds = (int) mxGetScalar(prhs[3]);
	int numColumnsShape = (int) mxGetScalar(prhs[4]);
	int numColumnsPred = (int) mxGetScalar(prhs[5]);
    int *colA = (int*) mxGetData(prhs[6]);
    double *colAFac = mxGetPr(prhs[7]);
    int *colB = (int*) mxGetData(prhs[8]);
    double *colBFac = mxGetPr(prhs[9]);
  	int *boundsPred = (int*) mxGetData(prhs[10]);

	/* intern variables and pointers */
	float* mu_a_b2 = NULL;
	float factor1,factor2,muMean;
	int i,j,k,idx,idxA,idxB;

  	/* 2-D matrix with [numBounds,numColumnsPred] */
	plhs[0] = mxCreateNumericMatrix(numRows,(numBounds-1)*numColumnsPred,mxSINGLE_CLASS,mxREAL);
    mu_a_b2 = (float *) mxGetPr(plhs[0]);

	/* negative entropy of q_c */
	#pragma omp parallel for private(k,i,idx,factor1,factor2,muMean,idxA,idxB)
	for (j=0; j < numColumnsPred; j++) {
		for (k=0; k < numBounds-1; k++) {
			idxA = colA[j]+k*numColumnsShape; idxB = colB[j]+k*numColumnsShape;
			muMean = (float) (colAFac[j]*mu[idxA+numColumnsShape] + colBFac[j]*mu[idxB+numColumnsShape]);
			factor1 = (float) (colAFac[j]*prec[idxA] + colBFac[j]*prec[idxB]);
			factor2 = (float) (colAFac[j]*mu[idxA]*prec[idxA] + colBFac[j]*mu[idxB]*prec[idxB]);
			idx = (k*numColumnsPred + j)*numRows;
			for (i=boundsPred[j*2]; i <= boundsPred[j*2+1]; i++) {
				mu_a_b2[idx + i] = muMean - (factor1*(i+1) - factor2);
			}
		}
	}
}
