#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

/* q_c.singleton = optQCMFC(condQB,prediction,Sigma_c,mu_c,c_c,mu_a_b,numColumnsPred,numColumnsShape,columnsPredShapeVec,columnsPredShapeFactorVec);
 * */

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

float expf_fast(float a) {
  union { float f; int x; } u;
  u.x = (int) (12102203 * a + 1064866805);
  return u.f;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	/* Input variables */
	double *condQB = mxGetPr(prhs[0]);
	double *prediction = mxGetPr(prhs[1]);
	double *mu_a_b = mxGetPr(prhs[2]);
	int numColumnsPred = (int) mxGetScalar(prhs[3]);
	int numColumnsShape = (int) mxGetScalar(prhs[4]);
	int *colA = (int*) mxGetData(prhs[5]);
	double *colAFac = mxGetPr(prhs[6]);
	int *colB = (int*) mxGetData(prhs[7]);
	double *colBFac = mxGetPr(prhs[8]);
	double *factorsPrec = mxGetPr(prhs[9]);
	double *hashTable = mxGetPr(prhs[10]);
	int *cmin = (int*) mxGetData(prhs[11]);
	int *cmax = (int*) mxGetData(prhs[12]);
	/* intern variables and pointers */
	double* q_c = NULL;
	double* boundaries = NULL;
	int i,j,k,i1,i2;
	int numRows = mxGetM(prhs[0]);
	int numBounds = mxGetN(prhs[0])/numColumnsPred;
	int alphaSize = numRows*numBounds*sizeof(double);
	double* alpha = malloc(alphaSize);
	double* beta = malloc(alphaSize); 

	double* c = malloc(numBounds*sizeof(double));
	double alphaTotal,q_c_total,tmpA,tmpB,valA,valB,factorA,factorB,cInv;
	double* preCalcA = malloc(numRows*sizeof(double));
	double* preCalcB = malloc(numRows*sizeof(double));

	int idxQC,idxPred,idx,idxA,idxB,idxC,idxANumRows,idxBNumRows,idxCond,idxBounds;
	int* A = malloc(numBounds*sizeof(int));
	int* B = malloc(numBounds*sizeof(int));

	/* determines when to use the hash table */
	/*int limit = 10;*/
	int limit2 = -30;
/*	int counter = 0;*/
	/* switch from matlab indexing to C indexing */
	for (i=0; i < numColumnsPred; i++) {
		colA[i] = colA[i]-1;
		colB[i] = colB[i]-1;
	}

	plhs[0] = mxCreateDoubleMatrix(1,numRows*numBounds*numColumnsPred,mxREAL);
	q_c = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(numBounds,numColumnsPred,mxREAL);
	boundaries = mxGetPr(plhs[1]);

	/* ****** start sum-product ******** */
	for (j=0; j < numColumnsPred; j++) {
/*	for (j=0; j < 1; j++) { */
		memset(alpha, 0, alphaSize);
		memset(beta, 0, alphaSize);

		idxPred = j;
		
		/* calculate limits of for-loops corresponding to transition matrices */
		for (k=0; k < numBounds; k++) {
			A[k] = cmin[idxPred + k*numColumnsPred]; 
			B[k] = cmax[idxPred + k*numColumnsPred];
			/*printf("%d, %d: %d, %d\n",j,k,A[k],B[k]);*/
		}
	
		alphaTotal = 0;
		/* shape index for condQB */
		idxA = colA[idxPred]*numRows; idxB = colB[idxPred]*numRows;
		/* pred index for prediction */
		idxC = idxPred*numRows*numBounds;
		for (i = A[0]; i <= B[0]; i++) {
			alpha[i] = condQB[j*numRows + i]*prediction[idxC + i];
			alphaTotal += alpha[i];
		}
		c[0] = alphaTotal; alphaTotal = 1/alphaTotal;

		/* normalize alpha */
		for (i=A[0]; i <= B[0]; i++) {
			alpha[i] *= alphaTotal;
		}

		/* make forward message passing over all boundaries */		
		/* for boundaries 2 to numBounds */
		for (k=1; k < numBounds; k++) {
/*		for(k=1; k < 0; k++) { */
			/* preCalc index for inner loop */
			idxA = colA[idxPred] + (k-1)*numColumnsShape; idxB = colB[idxPred] + (k-1)*numColumnsShape;
			factorA = -0.5*factorsPrec[idxA]; factorB = -0.5*factorsPrec[idxB];
			
			idxANumRows = idxA*numRows; idxBNumRows = idxB*numRows;
			idxCond = k*numColumnsPred*numRows + j*numRows;
			idx = numRows*k;

			alphaTotal = 0; 
			/* iterates over the columns of each transition matrix; corresponds to idxNonZeroA in matlab; determines the non-zero entries of the current alpha */
			#pragma omp parallel for private(tmpA,tmpB,valA,valB,i2) reduction(+:alphaTotal) 
			for (i1 = A[k]; i1 <= B[k]; i1++) {
				tmpA = tmpB = 0;
				/* iterates over the rows of transition matrices; corresponds to idxNonZeroB in matlab */
				/* upper triangular matrix --> ordering constraint on boundaries */
				for (i2 = A[k-1]; i2 <= min(i1,B[k-1]); i2++) {
					valA = (i1 + 1 - mu_a_b[idxANumRows + i2]);
					valA = valA*valA*factorA;
					valB = (i1 + 1 - mu_a_b[idxBNumRows + i2]);
					valB = valB*valB*factorB;
									
					if (valA > limit2) {tmpA += alpha[idx - numRows + i2]*hashTable[(int)(-valA*1000 + 0.5)];}
					if (valB > limit2) {tmpB += alpha[idx - numRows + i2]*hashTable[(int)(-valB*1000 + 0.5)];}
				}
				alpha[idx + i1] = prediction[idxC + idx + i1]*condQB[idxCond+i1]*(colAFac[idxPred]*tmpA + colBFac[idxPred]*tmpB);
				alphaTotal += alpha[idx + i1];
			}
			c[k] = alphaTotal; alphaTotal = 1/alphaTotal;

			/* normalize alpha */
			for (i = A[k]; i <= B[k]; i++) {
				alpha[idx + i] *= alphaTotal;
			}
		} /* end for over bounds k */

		/* init beta for the last node */
		idxQC = j*numBounds*numRows;
		idxBounds = (j+1)*numBounds - 1;
		boundaries[idxBounds] = 0;
		for (i=(numBounds-1)*numRows;i<numRows*numBounds;i++) {
			beta[i] = 1;
			q_c[idxQC + i] = alpha[i];
			boundaries[idxBounds] += alpha[i]*((i+1)-(numBounds-1)*numRows);
		}

		/* message backward */
		for (k=numBounds-2; k >= 0; k--) {
/*		for (k = 0; k < 0; k++) {*/
			idxCond = j*numRows + (k+1)*numColumnsPred*numRows;
			idxB = numRows*(k+1);	
			idxA = idxPred*numRows*numBounds + (k+1)*numRows;
			/* precalculate entries for inner loop over z_{n+1}, that are independent of z_n */
			for (i=A[k+1]; i <= B[k+1]; i++) {
				preCalcA[i] = beta[idxB + i]*prediction[idxA + i]*condQB[idxCond + i];
				/*preCalcB[i] = beta[idxB + i]*prediction[idxA + i]*condQB[idxCondA + i];*/
			}
			/* preCalc idx for inner loop */
			idxA = colA[idxPred] + k*numColumnsShape; idxB = colB[idxPred] + k*numColumnsShape;
			factorA = -0.5*factorsPrec[idxA]; factorB = -0.5*factorsPrec[idxB];
			
			idxANumRows = idxA*numRows; idxBNumRows = idxB*numRows; idx = numRows*k;

			/* the outer loop (over z_n) is constrained  by alpha (and therefor condQB), the inner loop over (z_{n+1}) by condQB */
			q_c_total = 0; cInv = 1/c[k+1];
			#pragma omp parallel for private(tmpA,tmpB,valA,valB,i2) reduction(+:q_c_total)
			for (i1 = A[k]; i1 <= B[k]; i1++) {
				tmpA = tmpB = 0;
				/* idxFinal */
				for (i2 = max(A[k+1],i1); i2 <= B[k+1]; i2++) {
					valA = factorA*(i2 + 1 - mu_a_b[idxANumRows + i1])*(i2 + 1 - mu_a_b[idxANumRows + i1]);
					valB = factorB*(i2 + 1 - mu_a_b[idxBNumRows + i1])*(i2 + 1 - mu_a_b[idxBNumRows + i1]);
			
				    if (valA > limit2) {tmpA += preCalcA[i2]*hashTable[(int)(-valA*1000 + 0.5)];}
                    if (valB > limit2) {tmpB += preCalcA[i2]*hashTable[(int)(-valB*1000 + 0.5)];}
				}
				beta[idx + i1] = (colAFac[idxPred]*tmpA + colBFac[idxPred]*tmpB)*cInv;
				q_c[idxQC + idx + i1] = alpha[idx + i1]*beta[idx + i1];
				q_c_total += q_c[idxQC + idx + i1];
			}
			idxBounds = j*numBounds + k;
			boundaries[idxBounds] = 0;
			/* convert to inverse */
			q_c_total = 1/q_c_total;
			/* normalize q_c distribution */
			for (i1 = A[k]; i1 <= B[k]; i1++) {
				q_c[idxQC + idx + i1] *= q_c_total;
				boundaries[idxBounds] += q_c[idxQC + idx + i1]*(i1+1);
			}
		}
	}
	free(alpha); free(beta); free(c); free(preCalcA); free(preCalcB); free(A); free(B);
}
