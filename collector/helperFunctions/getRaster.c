#include <mex.h>
#include <math.h>
#include <stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i,j,k;

	/* interpolation input with added first row zeros on last row ones */
	double *borders = mxGetPr(prhs[0]);
	/* dimension of the scan B0 */
	double* dim = mxGetPr(prhs[1]);

	double* Raster = NULL;
	int numBorders = mxGetM(prhs[0]);
	int idxStart,idxEnd;

	plhs[0] = mxCreateDoubleMatrix(dim[0],dim[1],mxREAL);
	/* output Raster matrix */
	Raster = mxGetPr(plhs[0]);

	for (j=0; j < (int)dim[1]; j++) {
		for (i=0;i < numBorders-1; i++) {
			idxStart = (borders[i + j*numBorders] < 0) ? 1 : borders[i + j*numBorders]+1; 
			/* idxStart = (int)fmax(borders[i + j*numBorders],(float) 0)+1; */
			idxEnd = (int)borders[i+1 + j*numBorders];
			if (idxStart<=idxEnd) {
				for (k=idxStart;k<=idxEnd;k++) {
					Raster[k-1 + j*(int)dim[0]] = i+1;
				}
			} 
		}
	} 
}
