#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

/* q_c.singleton = optQCMFC(condQB,prediction,Sigma_c,mu_c,c_c,mu_a_b,numColumnsPred,numColumnsShape,columnsPredShapeVec,columnsPredShapeFactorVec);
 * */

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

	clock_t begin, end;
	double time_elapsed_alpha,time_elapsed_beta;

	/* Counting variables */
	int numVolRegions = mxGetNumberOfElements(prhs[7]);
	int numRows = mxGetM(prhs[0]);

	/* determines when to use the hash table */
	int limit = -10;

	/* cumsum of the number of columns for previous regions */
	int numColsShapeCS[numVolRegions+1], ii;
	numColsShapeCS[0] = 0;
	for (ii=1;ii<numVolRegions+1;ii++) {
		numColsShapeCS[ii] = numColsShapeCS[ii-1] + (int)numColumnsShape[ii-1];
	}
	int numBounds = mxGetN(prhs[0])/numColsShapeCS[numVolRegions];
	plhs[0] = mxCreateDoubleMatrix(1,numRows*numBounds*numVolRegions*numColumnsPred,mxREAL);
	double *q_c= mxGetPr(plhs[0]);
	
	int i,j,k,i1,i2,vol;
	double alpha[numRows*numBounds];
/*	plhs[1] = mxCreateDoubleMatrix(1,numRows*numBounds,mxREAL);
	double *alpha = mxGetPr(plhs[1]); */
	double beta[numRows*numBounds]; 
/*	plhs[2] = mxCreateDoubleMatrix(1,numRows*numBounds,mxREAL);
	double *beta = mxGetPr(plhs[2]); */

	double c[numBounds];
	double alphaTotal,q_c_total,tmpA,tmpB,idxTmp,valA,valB,factorA,factorB;
	double preCalcA[numRows],preCalcB[numRows];

	int numColsPrevRegion,idxQC,idxShape,idxShapeWithout,idxPred,idx,idxA,idxB,idxC,idxANumRows,idxBNumRows,idxCondA,idxCondB;
	int idxCounter,idxCounterA,idxCounterB,numNonZeroIdx[numBounds];
	int idxNonZeroB[numRows], idxNonZeroAlpha[numBounds][numRows], numNonZeroIdxB;

	/* ****** start sum-product ******** */
	for (vol=0; vol < numVolRegions; vol++) {
	
		numColsPrevRegion = numColsShapeCS[vol]*numBounds;
		idxShape = numColsShapeCS[vol]*numBounds*numRows;
		idxShapeWithout = numColsShapeCS[vol]*(numBounds-1)*numRows;

		for (j=0; j < numColumnsPred; j++) {
			memset(alpha, 0, sizeof(alpha));
			memset(beta, 0, sizeof(beta));

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
					idxNonZeroAlpha[0][idxCounter] = i; idxCounter++;
				}
			}
			numNonZeroIdx[0] = idxCounter; idxCounter = 0;
			c[0] = alphaTotal;

			/* normalize */
			for (i=0; i < numNonZeroIdx[0]; i++) {
				alpha[idxNonZeroAlpha[0][i]] /= c[0];
			}

			/* for boundaries 2 to numBounds */
			for (k=1; k < numBounds; k++) {
				idxCounter = 0;
				idxA = idxShapeWithout +(k-1)*(int)numColumnsShape[vol]*numRows + ((int)colA[idxPred]-1)*numRows;
				idxB = idxShapeWithout +(k-1)*(int)numColumnsShape[vol]*numRows + ((int)colB[idxPred]-1)*numRows;
				idx = numRows*(k-1);
				/* precalculate entries for inner loop */
				for (idxCounter = 0; idxCounter < numNonZeroIdx[k-1]; idxCounter++) {
					i = idxNonZeroAlpha[k-1][idxCounter];
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
						idxNonZeroAlpha[k][idxCounterB] = i1; idxCounterB++;
						tmpA = tmpB = 0;
						/* iterates over the row of transition matrices; corresponds to idxNonZeroB in matlab */
						for (idxCounterA = 0; idxCounterA < numNonZeroIdx[k-1]; idxCounterA++) {
							i2 = idxNonZeroAlpha[k-1][idxCounterA];
							if (i2 <= i1 && preCalcA !=0 && preCalcB !=0) {
								valA = factorA*(i1 + 1 - mu_c[idxANumRows + i2])*(i1 + 1 - mu_c[idxANumRows + i2]);
								valB = factorB*(i1 + 1 - mu_c[idxBNumRows + i2])*(i1 + 1 - mu_c[idxBNumRows + i2]);
							
								/* printf("%d: %.2f\n",(int)(-valA*10 +0.5),valA); */
								tmpA += (valA < -limit) ? preCalcA[i2]*hashTable[(int)(-valA*100 +0.5)] : preCalcA[i2]*exp(valA);
								tmpB += (valB < -limit) ? preCalcB[i2]*hashTable[(int)(-valB*100 +0.5)] : preCalcB[i2]*exp(valB);
								
/*								tmpA += preCalcA[i2]*exp(-0.5*Sigma_inv_c[idxA]*pow(i1 + 1 - mu_c[idxANumRows + i2],2));
								tmpB += preCalcB[i2]*exp(-0.5*Sigma_inv_c[idxB]*pow(i1 + 1 - mu_c[idxBNumRows + i2],2)); */
							}
						}
						alpha[numRows*k + i1] = prediction[idxPred*numRows*numBounds + k*numRows + i1]*(colAFac[idxPred]*tmpA + colBFac[idxPred]*tmpB);
						alphaTotal += alpha[numRows*k + i1];
					}
				}

				c[k] = alphaTotal;
				/* time_elapsed_alpha += (double) ((double) (end - begin) / (double) CLOCKS_PER_SEC);*/

				/* copy idxNonZeroB onto idxNonZeroA and normalize alpha */
				for (i=0;i<idxCounterB; i++) {
					alpha[numRows*k + idxNonZeroAlpha[k][i]] /= c[k];
				}
				numNonZeroIdx[k] = idxCounterB;
 			}
		

			/* init beta */
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
					i1 = idxNonZeroAlpha[k][idxCounter];
					tmpA = tmpB = 0;
					/* idxFinal */
					for (idxCounterB = 0; idxCounterB < numNonZeroIdxB; idxCounterB++) {
						i2 = idxNonZeroB[idxCounterB];
						if (i2 >= i1) {
							valA = factorA*(i2 + 1 - mu_a_b[idxANumRows + i1])*(i2 + 1 - mu_a_b[idxANumRows + i1]);
							valB = factorB*(i2 + 1 - mu_a_b[idxBNumRows + i1])*(i2 + 1 - mu_a_b[idxBNumRows + i1]);
						
							/* printf("%d: %.2f\n",(int)(-valA*10 +0.5),valA); */
							tmpA += (valA < -limit) ? preCalcA[i2]*hashTable[(int)(-valA*100 +0.5)] : preCalcA[i2]*exp(valA);
							tmpB += (valB < -limit) ? preCalcB[i2]*hashTable[(int)(-valB*100 +0.5)] : preCalcB[i2]*exp(valB);
							/*tmpA += preCalcA[i2]*exp(-0.5*factorsPrec[idxA]*pow(i2 + 1 - mu_a_b[idxANumRows + i1] ,2));
							tmpB += preCalcB[i2]*exp(-0.5*factorsPrec[idxB]*pow(i2 + 1 - mu_a_b[idxBNumRows + i1] ,2)); */
						}
					}
					beta[idx + i1] = (colAFac[idxPred]*tmpA + colBFac[idxPred]*tmpB)/c[k+1];
					q_c[idxQC + idx + i1] = alpha[idx + i1]*beta[idx + i1];
					q_c_total += q_c[idxQC + idx + i1];
				}
				/* normalize q_c distribution */
				for (idxCounter = 0; idxCounter < numNonZeroIdx[k]; idxCounter++) {
					q_c[idxQC + idx + idxNonZeroAlpha[k][idxCounter]] /= q_c_total;
				}
			}
			/* time_elapsed_beta += (double) ((double) (end - begin) / (double) CLOCKS_PER_SEC); */
		}
	}
}
