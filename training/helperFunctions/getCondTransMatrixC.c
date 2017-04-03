#include <mex.h>
#include <math.h>
#include <stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	double* mu = mxGetPr(prhs[0]);
	double* prec = mxGetPr(prhs[1]);
	int numRows = (int) mxGetScalar(prhs[2]);
	double eps = (double) mxGetScalar(prhs[3]);
	int doTriu;
	if (nrhs > 4) {
		doTriu = (int) mxGetScalar(prhs[4]);
	} else {
		doTriu = 1;
	}
	/* only calculate probabilities up to this precision */
	double mu_a_b;
	int muFloor;
	double var_a_b = 1/prec[3]; /* for the conditional density, variance is given by the inverse precision matrix */
	double factor = pow(1/(2*3.1415926535897*var_a_b),0.5);
	double* i = NULL;
	double* j = NULL;
	double* s = NULL;
	int* numElements = NULL;
	int numBorders = mxGetM(prhs[0]);
	int boundA, boundB, startVal, stopVal;
	int numNotZero;
	int counter;
		
	/* calulate the number of elements larger than eps for each row */
	numNotZero = (int) ceil(sqrt(-log(eps*factor)*2*var_a_b))*2+1;
		
	plhs[0] = mxCreateDoubleMatrix(1,numNotZero*numRows,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1,numNotZero*numRows,mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1,numNotZero*numRows,mxREAL);
	plhs[3] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
	i = mxGetPr(plhs[0]);
	j = mxGetPr(plhs[1]);
	s = mxGetPr(plhs[2]);
	numElements = (int*) mxGetPr(plhs[3]);

	numNotZero = (numNotZero-1)/2;

	/* evaluate p(a|b) for all pairs of positions 1:numRows;
	 * Start at the rounded mean, and then run into both directions until we drop below precision "eps" */
	counter = 0;
	/* for positions of boundary k-1 (called here boundB): 1 till numRows-1 */
	for (boundB = 1; boundB <= numRows; boundB++) {
		mu_a_b = mu[1] - var_a_b*prec[1]*(boundB-mu[0]);
		muFloor = (int) mu_a_b;
		/* position of boundary k (called here boundA) */
		/* check for all possible values of boundB in this row: has to be at least boundA and can not be further away than numNotZero from muFloor */
		if (doTriu) {
			startVal = (muFloor-numNotZero < boundB) ? boundB : muFloor - numNotZero;
		} else {
			startVal = (muFloor-numNotZero < 1) ? 1 : muFloor - numNotZero;
		}
		stopVal = (muFloor+numNotZero > numRows) ? numRows : muFloor + numNotZero;
		for (boundA = startVal; boundA <= stopVal; boundA++) {
			i[counter] = boundB;
			j[counter] = boundA;
			s[counter] = factor*exp(-(boundA-mu_a_b)*(boundA-mu_a_b)*prec[3]*0.5);
			counter++;
		}
	}
	*numElements = counter;
}
