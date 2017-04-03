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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	/* Input variables */
	double *condQB = mxGetPr(prhs[0]);
	double *prediction = mxGetPr(prhs[1]);
	double *Sigma_inv_c = mxGetPr(prhs[2]);
	double *mu_c = mxGetPr(prhs[3]);
	double *c_c = mxGetPr(prhs[4]);
	double *mu_a_b = mxGetPr(prhs[5]);
	int numColumnsPred = (int) mxGetScalar(prhs[6]);
	int numColumnsShape = (int) mxGetScalar(prhs[7]);
	int *colA = (int*) mxGetData(prhs[8]);
	double *colAFac = mxGetPr(prhs[9]);
	int *colB = (int*) mxGetData(prhs[10]);
	double *colBFac = mxGetPr(prhs[11]);
	double *factorsPrec = mxGetPr(prhs[12]);
	double *hashTable = mxGetPr(prhs[13]);
	int *cmin = (int*) mxGetData(prhs[14]);
	int *cmax = (int*) mxGetData(prhs[15]);
	/* intern variables and pointers */
	double* q_c = NULL;

	int i,j,k,i1,i2;
	int numRows = mxGetM(prhs[0]);
	int numBounds = mxGetN(prhs[0])/numColumnsShape;
	int alphaSize = numRows*numBounds*sizeof(double);
	double* alpha = malloc(alphaSize);
	double* beta = malloc(alphaSize); 

	double* c = malloc(numBounds*sizeof(double));
	double alphaTotal,q_c_total,tmpA,tmpB,valA,valB,factorA,factorB,tmp;
	double* preCalcA = malloc(numRows*sizeof(double));
	double* preCalcB = malloc(numRows*sizeof(double));

	int idxQC,idxPred,idx,idxA,idxB,idxC,idxANumRows,idxBNumRows,idxCondA,idxCondB;
	int* A = malloc(numBounds*sizeof(int));
	int* B = malloc(numBounds*sizeof(int));

	/* determines when to use the hash table */
	int limit = 10;
	int limit2 = -30;
/*	int counter = 0;*/
	/* switch from matlab indexing to C indexing */
	for (i=0; i < numColumnsPred; i++) {
		colA[i] = colA[i]-1;
		colB[i] = colB[i]-1;
	}

	plhs[0] = mxCreateDoubleMatrix(1,numRows*numBounds*numColumnsPred,mxREAL);
	q_c = mxGetPr(plhs[0]);

	/* ****** start sum-product ******** */
	for (j=0; j < numColumnsPred; j++) {
	/*for (j=0; j < 1; j++) { */
		memset(alpha, 0, alphaSize);
		memset(beta, 0, alphaSize);

		idxPred = j;

		alphaTotal = 0;
		/* shape index for condQB */
		idxA = colA[idxPred]*numRows; idxB = colB[idxPred]*numRows;
		/* pred index for prediction */
		idxC = idxPred*numRows*numBounds;

		/* calculate limits of for-loops corresponding to transition matrices */
		for (k=0; k < numBounds; k++) {
			A[k] = min(cmin[colA[idxPred] + k*numColumnsShape],cmin[colB[idxPred] + k*numColumnsShape]);
			B[k] = min(cmax[colA[idxPred] + k*numColumnsShape],cmax[colB[idxPred] + k*numColumnsShape]);
			/*printf("%d, %d\n",A[k],B[k]);*/
		}
		for (i = A[0]; i <= B[0]; i++) {
			alpha[i] = (colAFac[idxPred]*condQB[idxA + i] + colBFac[idxPred]*condQB[idxB + i])*prediction[idxC + i];
			alphaTotal += alpha[i];
		}
		c[0] = alphaTotal;

		/* normalize alpha */
		for (i=A[0]; i <= B[0]; i++) {
			alpha[i] /= c[0];
		}

		/* make forward message passing over all boundaries */		
		/* for boundaries 2 to numBounds */
		for (k=1; k < numBounds; k++) {
			idxA = (k-1)*numColumnsShape*numRows + colA[idxPred]*numRows;
			idxB = (k-1)*numColumnsShape*numRows + colB[idxPred]*numRows;
			idx = numRows*(k-1);
			/* precalculate entries for inner loop */
			for (i = A[k-1]; i <= B[k-1]; i++) {
				preCalcA[i] = alpha[idx + i]*c_c[idxA + i];
				preCalcB[i] = alpha[idx + i]*c_c[idxB + i];
			}
			
			/* preCalc index for inner loop */
			idxA = colA[idxPred] + (k-1)*numColumnsShape; idxB = colB[idxPred] + (k-1)*numColumnsShape;
			idx = numRows*k;
			factorA = -0.5*Sigma_inv_c[idxA]; factorB = -0.5*Sigma_inv_c[idxB];
			
			idxANumRows = idxA*numRows; idxBNumRows = idxB*numRows;

			alphaTotal = 0; 
			/* iterates over the columns of each transition matrix; corresponds to idxNonZeroA in matlab; determines the non-zero entries of the current alpha */
			#pragma omp parallel for private(tmpA,tmpB,valA,valB,i2) reduction(+:alphaTotal)
			for (i1 = A[k]; i1 <= B[k]; i1++) {
				tmpA = tmpB = 0;
				/* iterates over the rows of transition matrices; corresponds to idxNonZeroB in matlab */
				/* upper triangular matrix --> ordering constraint on boundaries */
				for (i2 = A[k-1]; i2 <= min(i1,B[k-1]); i2++) {
					valA = factorA*(i1 + 1 - mu_c[idxANumRows + i2])*(i1 + 1 - mu_c[idxANumRows + i2]);
					valB = factorB*(i1 + 1 - mu_c[idxBNumRows + i2])*(i1 + 1 - mu_c[idxBNumRows + i2]);
				
					/*printf("%d, %d: %.2f, %d\n",i1,i2,valA,(int)(-valA*1000 + 0.5));*/
					/*if (valA > limit) { counter++; }; if (valB > limit) {counter++;}*/
					/*if (valA > limit2) {
						tmpA += (valA < limit) ? preCalcA[i2]*hashTable[(int)(-valA*1000 + 0.5)] : preCalcA[i2]*exp(valA);
					}
					if (valB > limit2) {
						tmpB += (valB < limit) ? preCalcB[i2]*hashTable[(int)(-valB*1000 + 0.5)] : preCalcB[i2]*exp(valB);
					}*/
					/*if (valA > limit2) {tmpA += preCalcA[i2]*exp(valA);}
					if (valB > limit2) {tmpB += preCalcB[i2]*exp(valB);}*/
					if (valA > limit2) {tmpA += preCalcA[i2]*hashTable[(int)(-valA*1000 + 0.5)];}
					if (valB > limit2) {tmpB += preCalcB[i2]*hashTable[(int)(-valB*1000 + 0.5)];}
				}
				alpha[idx + i1] = prediction[idxC + idx + i1]*(colAFac[idxPred]*tmpA + colBFac[idxPred]*tmpB);
				alphaTotal += alpha[idx + i1];
			}
			c[k] = alphaTotal;

			/* normalize alpha */
			for (i = A[k]; i <= B[k]; i++) {
				alpha[idx + i] /= c[k];
			}
		} /* end for over bounds k */

		/* init beta for the last node */
		idxQC = j*numBounds*numRows;
		for (i=(numBounds-1)*numRows;i<numRows*numBounds;i++) {
			beta[i] = 1;
			q_c[idxQC + i] = alpha[i];
		}

		/* message backward */
		for (k=numBounds-2; k >= 0; k--) {
			idxCondA = colA[idxPred]*numRows + (k+1)*numColumnsShape*numRows;
			idxCondB = colB[idxPred]*numRows + (k+1)*numColumnsShape*numRows;
			idxB = numRows*(k+1);	
			idxA = idxPred*numRows*numBounds + (k+1)*numRows;
			/* precalculate entries for inner loop over z_{n+1}, that are independent of z_n */
			for (i=A[k+1]; i <= B[k+1]; i++) {
				preCalcA[i] = beta[idxB + i]*prediction[idxA + i]*condQB[idxCondA + i];
				preCalcB[i] = beta[idxB + i]*prediction[idxA + i]*condQB[idxCondB + i];
			}
			/* preCalc idx for inner loop */
			idxA = colA[idxPred] + k*numColumnsShape; idxB = colB[idxPred] + k*numColumnsShape;
			factorA = -0.5*factorsPrec[idxA]; factorB = -0.5*factorsPrec[idxB];
			
			idxANumRows = idxA*numRows; idxBNumRows = idxB*numRows; idx = numRows*k;

			/* the outer loop (over z_n) is constrained  by alpha (and therefor condQB), the inner loop over (z_{n+1}) by condQB */
			q_c_total = 0;
			#pragma omp parallel for private(tmpA,tmpB,valA,valB,i2) reduction(+:q_c_total)
			for (i1 = A[k]; i1 <= B[k]; i1++) {
				tmpA = tmpB = 0;
				/* idxFinal */
				for (i2 = max(A[k+1],i1); i2 <= B[k+1]; i2++) {
					valA = factorA*(i2 + 1 - mu_a_b[idxANumRows + i1])*(i2 + 1 - mu_a_b[idxANumRows + i1]);
					valB = factorB*(i2 + 1 - mu_a_b[idxBNumRows + i1])*(i2 + 1 - mu_a_b[idxBNumRows + i1]);
			
					/*if (valA > limit) {	counter++;}; if (valB > limit) {counter++;}*/
/*					if (valA > limit2) {
						tmpA += (valA < limit) ? preCalcA[i2]*hashTable[(int)(-valA*1000 + 0.5)] : preCalcA[i2]*exp(valA);
					}
					if (valB > limit2) {
						tmpB += (valB < limit) ? preCalcB[i2]*hashTable[(int)(-valB*1000 + 0.5)] : preCalcB[i2]*exp(valB);
					}*/
				    /*if (valA > limit2) {tmpA += preCalcA[i2]*exp(valA);}
                    if (valB > limit2) {tmpB += preCalcB[i2]*exp(valB);}*/
				    if (valA > limit2) {tmpA += preCalcA[i2]*hashTable[(int)(-valA*1000 + 0.5)];}
                    if (valB > limit2) {tmpB += preCalcB[i2]*hashTable[(int)(-valB*1000 + 0.5)];}
				}
				beta[idx + i1] = (colAFac[idxPred]*tmpA + colBFac[idxPred]*tmpB)/c[k+1];
				q_c[idxQC + idx + i1] = alpha[idx + i1]*beta[idx + i1];
				q_c_total += q_c[idxQC + idx + i1];
			}
			/* normalize q_c distribution */
			for (i1 = A[k]; i1 <= B[k]; i1++) {
				q_c[idxQC + idx + i1] /= q_c_total;
			}
		}
	}
	free(alpha); free(beta); free(c); free(preCalcA); free(preCalcB); free(A); free(B);
	/* printf("counter: %d\n",counter); */
}
