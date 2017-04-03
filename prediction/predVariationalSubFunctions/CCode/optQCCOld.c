#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

/* q_c.singleton = optQCMFC(condQB,prediction,Sigma_c,mu_c,c_c,mu_a_b,numColumnsPred,numColumnsShape,columnsPredShapeVec,columnsPredShapeFactorVec);
 * */

/* sum the number of columns over all B-Scans */
int numberTotalColumns(double* numColumnsShape, int numVolRegions)
{
	int ii;
	int numColsTotal = 0;
	for (ii=0;ii<numVolRegions;ii++) {
		numColsTotal += (int)numColumnsShape[ii];
	}
	return numColsTotal;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	/* Input variables */
	double *condQB = mxGetPr(prhs[0]);
	double *prediction = mxGetPr(prhs[1]);
	double *Sigma_inv_c = mxGetPr(prhs[2]);
	double *mu_c = mxGetPr(prhs[3]);
	double *c_c = mxGetPr(prhs[4]);
	double *mu_a_b = mxGetPr(prhs[5]);
	int numColumnsPred = (int)(mxGetScalar(prhs[6]));
	double *numColumnsShape = mxGetPr(prhs[7]);
	double *colA = mxGetPr(prhs[8]);
	double *colAFac = mxGetPr(prhs[9]);
	double *colB = mxGetPr(prhs[10]);
	double *colBFac = mxGetPr(prhs[11]);
	double *factorsPrec = mxGetPr(prhs[12]);
	double *hashTable = mxGetPr(prhs[13]);

	/* intern variables and pointers */
	double* q_c = NULL;

	int i,j,k,i1,i2,vol;
	int numRows = mxGetM(prhs[0]);
	int numVolRegions = mxGetNumberOfElements(prhs[7]);
	int numBounds = mxGetN(prhs[0])/numberTotalColumns(numColumnsShape,numVolRegions);
	int alphaSize = numRows*numBounds*sizeof(double);
	double* alpha = malloc(alphaSize);
	double* beta = malloc(alphaSize); 

	double* c = malloc(numBounds*sizeof(double));
	double alphaTotal,q_c_total,tmpA,tmpB,valA,valB,factorA,factorB;
	double* preCalcA = malloc(numRows*sizeof(double));
	double* preCalcB = malloc(numRows*sizeof(double));

	int idxQC,idxShape,idxShapeWithout,idxPred,idx,idxA,idxB,idxC,idxANumRows,idxBNumRows,idxCondA,idxCondB;
	int idxCounter,idxCounterA,idxCounterB;
	int* numNonZeroIdx = malloc(numBounds*sizeof(int));
	int* idxNonZeroB = malloc(numRows*sizeof(int));
	/* int idxNonZeroAlpha[numBounds][numRows] */
   	int* idxNonZeroAlpha = malloc(numBounds*numRows*sizeof(int)); 
	int numNonZeroIdxB;

	/* determines when to use the hash table */
	int limit = -10;

	/* cumsum of the number of columns for previous regions */
	int* numColsShapeCS = malloc((numVolRegions+1)*sizeof(int)); 
	int ii;
	numColsShapeCS[0] = 0;
	for (ii=1;ii<numVolRegions+1;ii++) {
		numColsShapeCS[ii] = numColsShapeCS[ii-1] + (int)numColumnsShape[ii-1];
	}
	plhs[0] = mxCreateDoubleMatrix(1,numRows*numBounds*numVolRegions*numColumnsPred,mxREAL);
	q_c = mxGetPr(plhs[0]);
	

	/* ****** start sum-product ******** */
	for (vol=0; vol < numVolRegions; vol++) {
		idxShape = numColsShapeCS[vol]*numBounds*numRows;
		idxShapeWithout = numColsShapeCS[vol]*(numBounds-1)*numRows;

		for (j=0; j < numColumnsPred; j++) {
			memset(alpha, 0, alphaSize);
			memset(beta, 0, alphaSize);

			idxPred = numColumnsPred*vol + j;

			/* make forward message passing over all boundaries */		
			alphaTotal = 0;
			idxA = idxShape + ((int)colA[idxPred]-1)*numRows;
			idxB = idxShape + ((int)colB[idxPred]-1)*numRows;
			idxC = idxPred*numRows*numBounds;
			idxCounter = 0;
			for (i=0; i < numRows; i++) {
				if (condQB[idxA + i]!=0 || condQB[idxB + i] !=0) {
					alpha[i] = (colAFac[idxPred]*condQB[idxA + i] + colBFac[idxPred]*condQB[idxB + i])*prediction[idxC + i];
					alphaTotal += alpha[i];
					idxNonZeroAlpha[idxCounter] = i; idxCounter++;
				}
			}
			numNonZeroIdx[0] = idxCounter; idxCounter = 0;
			c[0] = alphaTotal;

			/* normalize */
			for (i=0; i < numNonZeroIdx[0]; i++) {
				alpha[idxNonZeroAlpha[i]] /= c[0];
			}

			/* for boundaries 2 to numBounds */
			for (k=1; k < numBounds; k++) {
				idxCounter = 0;
				idxA = idxShapeWithout +(k-1)*(int)numColumnsShape[vol]*numRows + ((int)colA[idxPred]-1)*numRows;
				idxB = idxShapeWithout +(k-1)*(int)numColumnsShape[vol]*numRows + ((int)colB[idxPred]-1)*numRows;
				idx = numRows*(k-1);
				/* precalculate entries for inner loop */
				for (idxCounter = 0; idxCounter < numNonZeroIdx[k-1]; idxCounter++) {
					i = idxNonZeroAlpha[(k-1)*numRows + idxCounter];
					preCalcA[i] = alpha[idx + i]*c_c[idxA + i];
					preCalcB[i] = alpha[idx + i]*c_c[idxB + i];
				}
				
				/* preCalc index for inner loop */
				idxA = idxShapeWithout/numRows + (int)colA[idxPred]-1 + (k-1)*(int)numColumnsShape[vol];
				idxB = idxShapeWithout/numRows + (int)colB[idxPred]-1 + (k-1)*(int)numColumnsShape[vol];
				factorA = -0.5*Sigma_inv_c[idxA]; factorB = -0.5*Sigma_inv_c[idxB];
			    	
				idxANumRows = idxA*numRows; idxBNumRows = idxB*numRows;

				/* preCalc idx for outer loop */
				idxCondA = idxShape + ((int)colA[idxPred]-1)*numRows + k*(int)numColumnsShape[vol]*numRows;
				idxCondB = idxShape + ((int)colB[idxPred]-1)*numRows + k*(int)numColumnsShape[vol]*numRows;

				idxCounterB = 0; alphaTotal = 0;
				/* iterates over the columns of each transition matrix; corresponds to idxNonZeroA in matlab; determines the non-zero entries of the current alpha */
				for (i1=0; i1 < numRows; i1++) {
					if (condQB[idxCondA + i1] != 0 || condQB[idxCondB + i1] != 0) {
						idxNonZeroAlpha[k*numRows + idxCounterB] = i1; idxCounterB++;
						tmpA = tmpB = 0;
						/* iterates over the rows of transition matrices; corresponds to idxNonZeroB in matlab */
						for (idxCounterA = 0; idxCounterA < numNonZeroIdx[k-1]; idxCounterA++) {
							i2 = idxNonZeroAlpha[(k-1)*numRows + idxCounterA];
							if (i2 <= i1 && preCalcA !=0 && preCalcB !=0) {
								valA = factorA*(i1 + 1 - mu_c[idxANumRows + i2])*(i1 + 1 - mu_c[idxANumRows + i2]);
								valB = factorB*(i1 + 1 - mu_c[idxBNumRows + i2])*(i1 + 1 - mu_c[idxBNumRows + i2]);
							
								/* printf("%d: %.2f\n",(int)(-valA*10 +0.5),valA); */
								if (valA > -300) {
									tmpA += (valA < -limit) ? preCalcA[i2]*hashTable[(int)(-valA*1000 +0.5)] : preCalcA[i2]*exp(valA);
								}
								if (valB > -300) {
									tmpB += (valB < -limit) ? preCalcB[i2]*hashTable[(int)(-valB*1000 +0.5)] : preCalcB[i2]*exp(valB);
								}
							}
						}
						alpha[numRows*k + i1] = prediction[idxPred*numRows*numBounds + k*numRows + i1]*(colAFac[idxPred]*tmpA + colBFac[idxPred]*tmpB);
						alphaTotal += alpha[numRows*k + i1];
					}
				}
				c[k] = alphaTotal;

				/* copy idxNonZeroB onto idxNonZeroA and normalize alpha */
				for (i=0;i<idxCounterB; i++) {
					alpha[numRows*k + idxNonZeroAlpha[k*numRows + i]] /= c[k];
				}
				numNonZeroIdx[k] = idxCounterB;
 			}

			/* init beta for the last node */
			idxQC = (vol*(numColumnsPred*numBounds) + j*numBounds)*numRows;
			for (i=(numBounds-1)*numRows;i<numRows*numBounds;i++) {
				beta[i] = 1;
				q_c[idxQC + i] = alpha[i];
			}

			/* message backward */
			for (k=numBounds-2; k >= 0; k--) {
				idxCondA = idxShape + ((int)colA[idxPred]-1)*numRows + (k+1)*(int)numColumnsShape[vol]*numRows;
				idxCondB = idxShape + ((int)colB[idxPred]-1)*numRows + (k+1)*(int)numColumnsShape[vol]*numRows;
				idxB = numRows*(k+1);	
				idxA = idxPred*numRows*numBounds + (k+1)*numRows;
				idxCounter = 0;
				/* precalculate entries for inner loop*/
				for (i=0; i < numRows; i++) {
					if (beta[idxB+i]!=0 && (condQB[idxCondA + i] != 0 || condQB[idxCondB + i] != 0)) {
						idxNonZeroB[idxCounter] = i; idxCounter++;
						preCalcA[i] = beta[idxB + i]*prediction[idxA + i]*condQB[idxCondA + i];
						preCalcB[i] = beta[idxB + i]*prediction[idxA + i]*condQB[idxCondB + i];
					}
				}
				numNonZeroIdxB = idxCounter;
				/* preCalc idx for inner loop */
				idxA = idxShapeWithout/numRows + (int)colA[idxPred]-1 + k*(int)numColumnsShape[vol];
			    idxB = idxShapeWithout/numRows + (int)colB[idxPred]-1 + k*(int)numColumnsShape[vol];
				factorA = -0.5*factorsPrec[idxA]; factorB = -0.5*factorsPrec[idxB];
				idxANumRows = idxA*numRows; idxBNumRows = idxB*numRows; idx = numRows*k;

				/* idxB */
				q_c_total = 0;
				for (idxCounter = 0; idxCounter < numNonZeroIdx[k]; idxCounter++) {
					i1 = idxNonZeroAlpha[k*numRows + idxCounter];
					tmpA = tmpB = 0;
					/* idxFinal */
					for (idxCounterB = 0; idxCounterB < numNonZeroIdxB; idxCounterB++) {
						i2 = idxNonZeroB[idxCounterB];
						if (i2 >= i1) {
							valA = factorA*(i2 + 1 - mu_a_b[idxANumRows + i1])*(i2 + 1 - mu_a_b[idxANumRows + i1]);
							valB = factorB*(i2 + 1 - mu_a_b[idxBNumRows + i1])*(i2 + 1 - mu_a_b[idxBNumRows + i1]);
					
							if (valA > -300) {
								tmpA += (valA < -limit) ? preCalcA[i2]*hashTable[(int)(-valA*1000 + 0.5)] : preCalcA[i2]*exp(valA);
							}
							if (valB > -300) {
								tmpB += (valB < -limit) ? preCalcB[i2]*hashTable[(int)(-valB*1000 + 0.5)] : preCalcB[i2]*exp(valB);
							}
						}
					}
					beta[idx + i1] = (colAFac[idxPred]*tmpA + colBFac[idxPred]*tmpB)/c[k+1];
					q_c[idxQC + idx + i1] = alpha[idx + i1]*beta[idx + i1];
					q_c_total += q_c[idxQC + idx + i1];
				}
				/* normalize q_c distribution */
				for (idxCounter = 0; idxCounter < numNonZeroIdx[k]; idxCounter++) {
					q_c[idxQC + idx + idxNonZeroAlpha[k*numRows + idxCounter]] /= q_c_total;
				}
			}
		}
	}
	free(alpha); free(beta); free(c); free(preCalcA); free(preCalcB); free(numNonZeroIdx); free(idxNonZeroB); free(idxNonZeroAlpha); free(numColsShapeCS);
}
