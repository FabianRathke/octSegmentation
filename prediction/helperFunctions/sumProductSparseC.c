#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	double* prediction = mxGetPr(prhs[0]); /* containts evaluated appearance models for all pixel in that BScan and all boundary models */
	double* mu = mxGetPr(prhs[1]);
	double* WML = mxGetPr(prhs[2]);
	double sigmaML = (double) mxGetScalar(prhs[3]);
	int* idxA = (int*) mxGetData(prhs[4]);
  	double *hashTable = mxGetPr(prhs[5]);

/*	int* idxRows = (int*) mxGetPr(prhs[4]); */
	int M = mxGetN(prhs[2]); /* number of eigenmodes of the shape prior */
	int N = mxGetM(prhs[2]); /* number of rows in WML */
	int numColumns = mxGetNumberOfElements(prhs[4]); /* for each column in the shape prior model there is one index */
	const int* dims = mxGetDimensions(prhs[0]);
	int numBounds = dims[1]; /* number of columns of pObs */
	int numRows = mxGetM(prhs[0]); /* number of rows of pObs */
 	int numColumnsShape = (int) (double)N/(double)numBounds; /* for each column in the shape prior model there is one index */


	/* only calculate pairiwse probabilities up to this precision */
	double eps = pow(10,-30);
	double* gamma = NULL;
	/* return the probabilities for this B-Scan column */
	plhs[0] = mxCreateDoubleMatrix(1,numColumns*numBounds*numRows,mxREAL);
	gamma = mxGetPr(plhs[0]);

	#pragma omp parallel
	{	
		double var_a_b, factor, mu_a, mu_b, prec_a_b, prec_a_a, prec_b_b, evalDens, aTilde;
		double variance, varInv;
		double* prec = malloc(2*(numBounds-1)*sizeof(double));
		double var1, var2, var3;

		double* alpha = malloc(numBounds*numRows*sizeof(double));
		double* beta = malloc(numBounds*numRows*sizeof(double));
		double* c = malloc(numBounds*sizeof(double));
		double cInv, cTmp;
		double* prodTmp = malloc(numRows*sizeof(double));
		int boundA, boundB;
		int numNotZero;
		int i,k,idx,column,predOffset;
		int muFloor, startVal, stopVal;
		double mu_b_a, mu_a_b, tmpVal, prec_a_aaTilde;
		#pragma omp for	
		for (column=0; column < numColumns; column++) {
			memset(alpha,0,numRows*numBounds*sizeof(double));
			memset(beta,0,numRows*numBounds*sizeof(double));
			predOffset = numBounds*numRows*column; /* shifts the pointer inside the prediction matrix to the next BScan column */
			/* calculate for the first boundary; calculate the 1-d shape marginal p(a) for column j on the fly  */
			/* calculate the variance for the shape prior density: \Sigma = WW^T + sigma^2I */
			variance = 0;
			for (i=0; i < M; i++) {
				variance += WML[idxA[column] + i*N]*WML[idxA[column] + i*N];
			}
			
			variance += sigmaML; varInv = -0.5/variance;
			factor = 1/sqrt(2*3.1415926535897*variance);

			cTmp = 0;
			for (i=0; i < numRows; i++) {
				alpha[i] = factor*exp(varInv*(i+1 - mu[idxA[column]])*(i+1-mu[idxA[column]]))*prediction[i + predOffset];
				cTmp += alpha[i];
			}

			c[0] = cTmp; cInv = 1/cTmp;
			for (i=0; i < numRows; i++) {
				alpha[i] = alpha[i]*cInv;
			}

			/* calculate the precision matrices required for conditional densities p(a|b) */
			for (k=0; k < numBounds-1; k++) {
				var1 = 0; var2 = 0; var3 = 0;
				for (i=0; i < M; i++) {
					var1 += WML[idxA[column] + numColumnsShape*k + i*N]*WML[idxA[column] + numColumnsShape*k + i*N];
					var2 += WML[idxA[column] + numColumnsShape*k + i*N]*WML[idxA[column] + numColumnsShape*(k+1) + i*N];
					var3 += WML[idxA[column] + numColumnsShape*(k+1) + i*N]*WML[idxA[column] + numColumnsShape*(k+1) + i*N];
				}
				var1 += sigmaML;
				var3 += sigmaML;

				factor = 1/(var1*var3 - var2*var2);
				prec[0 + k*2] = factor*(-var2);
				prec[1 + k*2] = factor*var1;
			}

			/* calculate the remaining boundaries */
			for (k=1; k < numBounds; k++) {
				/* calculate parameters of distribution p(b|a); b is the boundary k-1; a is boundary k */
				var_a_b = 1/prec[1 + (k-1)*2]; /* for the conditional density, variance is given by the inverse precision matrix */
				factor = 1/sqrt(2*3.1415926535897*var_a_b);
				/* calulate the number of elements larger than eps for each row */
				cTmp = 0;
				mu_b = mu[idxA[column] + numColumnsShape*(k-1)]; mu_a = mu[idxA[column] + numColumnsShape*k]; prec_a_b = prec[0 + (k-1)*2]; prec_a_a = prec[1 + (k-1)*2]; aTilde = prec_a_b/prec_a_a;
				prec_a_aaTilde = 0.5*(1/prec_a_a)*prec_a_b*prec_a_b;
	 			numNotZero = (int) ceil(abs(sqrt(-log(eps*factor)*2*var_a_b)/aTilde));
				/* value of boundary k */
				for (boundA = 1; boundA <= numRows; boundA++) {
					mu_b_a = ((mu_a - boundA)/aTilde + mu_b);
					muFloor = (int) mu_b_a;
					/* the position of boundary k-1 can lie between 1 and boundary k and is constrained to lie within 2*numNotZero + 1 around its mean mu_b_a */
					startVal = (muFloor-numNotZero < 1) ? 1 : muFloor - numNotZero;
					stopVal = (muFloor+numNotZero > boundA) ? boundA : muFloor + numNotZero;
					tmpVal = 0;
					/* value of boundary k-1 */
					for (boundB = startVal; boundB <= stopVal; boundB++) {
						/* evaluate the conditional density; but only the part inside the exp function */
						evalDens = -(boundB-mu_b_a)*(boundB-mu_b_a)*prec_a_aaTilde;
						if (evalDens < -10) {
							evalDens = hashTable[(int)(-evalDens*1000 + 0.5)];
						} else {
							evalDens = exp(evalDens);
						}
						tmpVal += factor*evalDens*alpha[boundB-1+numRows*(k-1)];
					}
					alpha[boundA-1+numRows*k] = tmpVal*prediction[boundA-1+numRows*k + predOffset];
					cTmp += alpha[boundA-1+numRows*k];
				}
				/* normalize alpha */
				c[k] = cTmp; cInv = 1/cTmp;
				for (i=0; i < numRows; i++) {
					alpha[i+numRows*k] = alpha[i+numRows*k]*cInv;
				}
			}

			/* initialize beta */
			for (i=0; i < numRows; i++) {
				idx = (numBounds-1)*numRows+i;
				beta[idx] = 1;
				gamma[idx + predOffset] = alpha[idx]*beta[idx];
			}

		   for (k=numBounds-1; k > 0; k--) {
				/* precalculate the product of pObs*beta(z_{n+1}) */
				cInv = 1/c[k];
				for (i=0; i < numRows; i++) {
					prodTmp[i] = cInv*prediction[numRows*k+i + predOffset]*beta[numRows*k+i];
				}

				var_a_b = 1/prec[1 + (k-1)*2]; /* for the conditional density, variance is given by the inverse precision matrix */
				factor = 1/sqrt(2*3.1415926535897*var_a_b);
				/* calulate the number of elements larger than eps for each row */
				numNotZero = (int) ceil(sqrt(-log(eps*factor)*2*var_a_b));
				cTmp = 0;
				mu_b = mu[idxA[column] + numColumnsShape*(k-1)]; mu_a = mu[idxA[column] + numColumnsShape*k]; prec_a_b = prec[0 + (k-1)*2]; prec_a_a = prec[1 + (k-1)*2]*0.5;
				for (boundB = 1; boundB <= numRows; boundB++) {
					mu_a_b = mu_a - var_a_b*prec_a_b*(boundB-mu_b);
					muFloor = (int) mu_a_b;
					/* position of boundary k (called here boundA) */
					/* check for all possible values of boundB in this row: has to be at least boundA, can be at most numRows and we limit it to be not further away from muFloor than numNotZero */
					startVal = (muFloor-numNotZero < boundB) ? boundB : muFloor - numNotZero;
					stopVal = (muFloor+numNotZero > numRows) ? numRows : muFloor + numNotZero;
					tmpVal = 0;
					for (boundA = startVal; boundA <= stopVal; boundA++) {
						evalDens = -(boundA-mu_a_b)*(boundA-mu_a_b)*prec_a_a;
						evalDens = hashTable[(int)(-evalDens*1000 + 0.5)];
				/*		evalDens = exp(evalDens); */
						tmpVal += factor*evalDens*prodTmp[boundA-1];
					}
					idx = boundB-1+numRows*(k-1);
					beta[idx] = tmpVal;
					gamma[idx + predOffset] = beta[idx]*alpha[idx];
				}
			}
		}
		free(alpha); free(beta); free(prodTmp); free(prec); free(c);
	}
}	
