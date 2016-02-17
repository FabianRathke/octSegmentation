#include <math.h>
#include "matrix.h" /* mx */
#include "mex.h"

/* Regularized covariance estimation, "Sparse inverse covariance estimation
 * with the graphical lasso", Friedman et al., 2008
 *
 * Matlab interface:
 * function [W P] = cglasso2(S, lambda, penDiag, zero_elements_P, zero_elements_W)
 *
 * INPUT
 * ***************
 * 	S 			empirical covariance matrix
 *	lambda			regularization strength parameter
 *	penDiag			if true, then add lambda to the diag of S
 *	zero_elements_P	zero/one matrix of size(S) that forces entries of P
 *                       to become zero (if zero_elements(i,j)==1 --> P(i,j)=0)
 *	zero_elements_W	zero/one matrix of size(S) that forces entries of W
 *                       to become zero (if zero_elements(i,j)==1 --> W(i,j)=0)
 *
 * OUTPUT
 * **************
 *	W				estimated covariance matrix
 *	P				estimated precision matrix
 *
 */
#define DEBUG_PRINT(name) mexPrintf(#name " = %g\n", name);
#define DEBUG_INT(name) mexPrintf(#name " = %d\n", name);

const unsigned undefined_u = 12345;
const double undefined_d = 12345.6789;

/* p == matrix dimension */
#define a_(X, j, k, p) ((X)[(j) * (p) + (k)])
#define b_(X, j, k, p) ((X)[(j) * ((p) - 1) + (k)])

double mean_abs(double* beta, unsigned p) {
	double ma = 0;
	unsigned k = 0;
	for (k = 0; k < p - 1; ++k)
		ma += fabs(beta[k]);
	return ma / (p - 1); /* the other day, drop p - 1 from here and from "diff" */
}

void clean_zeros(double* zero_m, unsigned p) {
	/* set_to_zero = repmat(1:p,p,1); */
	/* if nargin > 3 */
	/* 	if ~isequal (zero_elements,zero_elements') */
	/* 		warning('Enforce symmetry of zero elements matrix'); */
	/* 		zero_elements = zero_elements + zero_elements'; */
	/* 		zero_elements(zero_elements>0) = 1; */
	/* 	end */
	/* 	set_to_zero(logical(zero_elements)) = 0; */
	/* end */
	unsigned j = undefined_u;
	unsigned k = undefined_u;
	unsigned diag_warning = 0;
	unsigned sparse_warning = 0;
	for (j = 0; j < p; ++j) {
		/* enforce diagonal == 0, although unused in theory. */
		double* z = &a_(zero_m, j, j, p);
		if (*z != 0) {
			diag_warning = 1;
			*z = 0;
		}
		for (k = 0; k < p; ++k) {
			double* z1 = &a_(zero_m, j, k, p);
			double* z2 = &a_(zero_m, k, j, p);
			if (*z1 != *z2)
				sparse_warning = 1;
			*z1 = *z2 = !!*z1 || !!*z2;
		}
	}
	if (diag_warning) 
		mexWarnMsgIdAndTxt("cglasso:warning",
			"Reset diagonal entry of zero elements matrix");

	if (sparse_warning)
		mexWarnMsgIdAndTxt("cglasso:warning",
			"Enforcing symmetry of zero elements matrix");
}

void copy_zeros(double* zero_ws, double* W, unsigned p) {
	double* W_end = W + p * p;
	for (; W != W_end; ++W)
		if (*zero_ws++)
			*W = 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	unsigned p = undefined_u; /* matrix dimension  */
	unsigned i = undefined_u; /* index */
	unsigned j = undefined_u; /* index */
	unsigned k = undefined_u; /* index */
	const unsigned ret_count = 2;
	const unsigned arg_count = 5;
	const unsigned arg_is_matrix[] = { 1, 0, 0, 1, 1 };
	const unsigned arg_is_double[] = { 1, 1, 0, 1, 1 };
	/* various index/etc. variables */
	unsigned t = undefined_u;
	const mxArray* S_array = NULL;
	const mxArray* zero_array_P = NULL;
	const mxArray* zero_array_W = NULL;
	mxArray* W_array = NULL;
	mxArray* beta_storage_array = NULL;
	mxArray* w_12_array = NULL;
	mxArray* P_array = NULL;
	mxArray* zero_betas_array = NULL;
	mxArray* zero_ws_array = NULL;
	const double* S = NULL;
	double* zero_betas = NULL;
	double* zero_ws = NULL;
	double* W = NULL;
	double* beta_storage = NULL;
	double* w_12 = NULL;
	double* P = NULL;
	double lambda = undefined_d;
	double penDiag = undefined_d;
	unsigned* skips = NULL;
	unsigned* skip = NULL;
	unsigned* skip2 = NULL;
	double inf[2] = {1, 0}; /* work around "/ 0" compiler warnings */
	double mean_abs_change = undefined_d;
	double mean_abs_change_norm = undefined_d;
	double threshold = undefined_d;
	const double* ts = NULL;
	double threshold_sub = undefined_d;
	double abs_change = undefined_d;
	double* beta = NULL;
	double diff = undefined_d;
	double old_beta = undefined_d;
	double new_beta = undefined_d;
	unsigned sparse_P = undefined_u;
	unsigned sparse_W = undefined_u;
	int iter;
	double x;
	int sign;

	/* first, check if arguments adhere to the input interface */
	if (nrhs != arg_count)
		mexErrMsgIdAndTxt("cglasso:syntax", "%d in place of %d "
				"arguments of cglasso2.", nrhs, arg_count);
	if (nlhs != ret_count)
		mexErrMsgIdAndTxt("cglasso:syntax", "%d in place of %d "
				"return values to cglasso2.", nlhs, ret_count);
	for (p = 0, t = 0; t != arg_count; ++t) {
		if (arg_is_matrix[t]) {
			/* test for quadratic matrix */
			const mxArray* x_array = prhs[t];
			unsigned x_size = mxGetN(x_array);
			if (x_size != mxGetM(x_array))
				mexErrMsgIdAndTxt("cglasso:syntax",
					"argument number %d of cglasso must be"
					" a quadratic matrix", t + 1);
			/* [p p] = size(S); */
			if (!p)
				p = x_size;
			/* other matrices must be either empty or
							of the same size p */
			else if (x_size && (x_size != p))
				mexErrMsgIdAndTxt("cglasso:syntax",
					"argument number %d of cglasso must be"
					" of dimension %d", t + 1, p);
		}
	}
	for (t = 0; t != arg_count; ++t)
		if (arg_is_double[t] && !mxIsDouble(prhs[t]))
			mexErrMsgIdAndTxt("cglasso:syntax", "argument number "
					"%d of cglasso must be of Matlab class"
					" 'double'", t + 1);

	S_array = prhs[0];
	S = mxGetPr(S_array);
	lambda = mxGetScalar(prhs[1]);
	penDiag = mxGetScalar(prhs[2]); /* automatically converts to double. */
	zero_array_P = prhs[3];
	sparse_P = (mxGetN(zero_array_P) > 1);
	zero_array_W = prhs[4];
	sparse_W = (mxGetN(zero_array_W) > 1);

	if (lambda < 0)
		mexErrMsgIdAndTxt("cglasso:syntax", "negative lambda");

/* if lambda == 0  */
/* 	W = S;  */
/* 	P = inv(W); */
/* 	return; */
/* end */
	if (lambda == 0)
		mexErrMsgIdAndTxt("cglasso:unimplemented", "lambda == 0");

	/* % initialize variables */

	/* allocate matrices W, P; W = S */
	W_array = mxDuplicateArray(S_array);
	W = mxGetPr(W_array);
	P_array = mxCreateDoubleMatrix(p, p, mxREAL);
	P = mxGetPr(P_array);

	if (sparse_P) {
		/* copy zero array to allow for modifications */
		/* (should be bits, not doubles.. */
		zero_betas_array = mxDuplicateArray(zero_array_P);
		zero_betas = mxGetPr(zero_betas_array);
	}
	if (sparse_W) {
		zero_ws_array = mxDuplicateArray(zero_array_W);
		zero_ws = mxGetPr(zero_ws_array);
	}
	
	if (penDiag)
		for (i = 0; i < p; ++i)
			a_(W, i, i, p) += lambda; /* W = S + diag(ones(1,p)*lambda); */

	/* beta_storage = zeros(p,p-1); */
	/* allocate matrix beta_storage */
	beta_storage_array = mxCreateDoubleMatrix(p, p - 1, mxREAL);
	beta_storage = mxGetPr(beta_storage_array);
	w_12_array = mxCreateDoubleMatrix(p - 1, 1, mxREAL);
	w_12 = mxGetPr(w_12_array);

	/* % for fast indexing */
	/* % *: idx == off-diagonal(p, p) */
	/* idx = logical(ones(p,p) - diag(ones(1,p))); */
	/* % *: idx2 == off-diagonal(p - 1, p - 1) */
	/* idx2 = logical(ones(p-1,p-1) - diag(ones(1,p-1))); */

	skips = (unsigned *) mxMalloc((sizeof(unsigned)) * (p - 1) * p);
	skip = skips;
	for (j = 0; j < p; ++j)
		for (k = 0; k < p; ++k)
			if (k != j)
				*skip++ = k;
	
	/* % stoping criterion for the outer loop */
	mean_abs_change =  inf[0] / inf[1]; /* mean_abs_change = inf; */
	mean_abs_change_norm = 1 / ((((double) p) - 1) * p)
	; /* .../(p^2-p); */

	/* threshold = 0.001 * mean(mean(abs(S-diag(diag(S))))); */
	/* sum(abs(off-diagonal)) / p^2 */
	threshold = 0;
	ts = S;
	for (k = 0; k < p - 1; ++k) {
		++ts;
		for (j = 0; j < p; ++j)
			threshold += fabs(*ts++);
	}
	threshold = 0.001 * threshold / p / p;
/* DEBUG_PRINT(threshold); */
	threshold_sub = 0.001;

	/* % check whether some elements of the precision matrix are supposed to be zero */
/*DEBUG_INT(sparse_P)
DEBUG_INT(sparse_W)*/
	if (sparse_P)
		clean_zeros(zero_betas, p);
	/* dito for W */
	if (sparse_W) {
		clean_zeros(zero_ws, p);
		copy_zeros(zero_ws, W, p);
	}

	/* % repeat block coordinate descent as long as the mean change of elements of W */
	/* % after one iteration is bigger than the mean abs value of entries of S  */
	/* %								times 0.001 */

	iter = 1;
	while (mean_abs_change > threshold && iter < 50) {
	iter++;

/* mexPrintf("mean_abs_change = %g\n", mean_abs_change); */ 
/* 		% pick a random sequence of columns */
/* 		randidx = my_randperm(p); */

		/* % mean absolute change of elements in W during one iteration, used */
		/* %						 as stopping criteria */
		abs_change = 0;
	 
		/* % sweep over all columns (block coordinate descent) */
		skip = skips; /* indexed by j */
		for (j = 0; j < p; ++j, skip += (p - 1)) { /* for j = randidx */
			/* % extract needed parts of blockmatrices */
			/* W_11 = W(idx(j,:),idx(j,:)); // remove row & column j */
			/* s_12 = S(j,idx(j,:));	// row j minus entry j */
			/* W_11[x, y] == a_(W, skip[x], skip[y]) */
			/* s_12[y] == a_(S, j, skip[y]) */
 
			/* % initialize beta */
			/* beta = beta_storage(j,:); */
			beta = beta_storage + j * (p - 1);
			diff = inf[0] / inf[1]; /* diff = inf; */
 
/* 			% entries that are supposed to be zero selected by the user; */
/* 			% compensate for missing diagonal entries by adding 1 */
/* 			set_to_zero_tmp = [find(set_to_zero(j,1:j-1)==0) ... */
/* 						find(set_to_zero(j,j+1:end)==0)+1]; */
/* 			% select entries of beta to iterate over */
 
			/* % solve the lasso subproblem by coordinate descent */
			while (diff > mean_abs(beta, p) * threshold_sub)
						/* mean(abs(beta)*threshold_sub) */
			{
/* mexPrintf("mean_abs(beta, p) = %g\n", mean_abs(beta, p)); */
/* 				beta_old = beta; */
				diff = 0; /* cumulate diff in-place */
/* 				% filter out elements that are defined to be zero */
/* 				%					by the user */
/* 				idx_sub = my_randperm(p-1); */
/* 				%idx_sub = setdiff(idx_sub,set_to_zero_tmp);  */
/* 				% ^*: don't actually implement the set-to-zero feature.. */

				skip2 = skips; /* indexed by k */
				for (k = 0; k < p - 1; ++k, skip2 += (p - 1)) { /* for k = idx_sub */
					/* % calc components needed for soft treshholding */
					/* x = s_12(k) ... */
					/*    - W_11(k,idx2(k,:))*beta(idx2(k,:))'; */
					if (sparse_P && a_(zero_betas, j, skip[k], p))
						continue;
					x = a_(S, j, skip[k], p);
					for (i = 0; i < p - 2; ++i) {
						x -= a_(W, skip[k], skip[skip2[i]], p)
							* beta[skip2[i]];
					}
					/* % make soft thresholding */
					/* beta(k) = sign(x)*max(abs(x)-lambda,0) ... */
					/* 				/ W_11(k,k); */
/*					oldCode: new_beta = (1 - 2 * !!signbit(x))
						* fmax(fabs(x) - lambda, 0)
						/ a_(W, skip[k], skip[k], p); */
					sign = (x>0) ? 1 : -1; 
					new_beta = (fabs(x)-lambda < 0) ? 0 : sign*(fabs(x)-lambda)/a_(W, skip[k], skip[k], p) ; 
					old_beta = beta[k];
					beta[k] = new_beta;
					/* % compute the difference between the old and the new */
					/* %						beta */
					/* diff = mean(abs(beta_old - beta)); */
					diff += fabs(old_beta - new_beta);
				}
				/* diff = mean(abs(beta_old - beta)); */
				diff /= (p - 1);
			}
			/* beta_storage(j,:) = beta; */
			/* w_12 = W_11*beta'; */
			/* naive matrix multiplication in C, hmm, hmm... */
			for (i = 0; i < p - 1; ++i) {
				w_12[i] = 0;
				for (k = 0; k < p - 1; ++k) {
					w_12[i] += a_(W, skip[i], skip[k], p) * beta[k];
				}
			}
			/* abs_change = abs_change + sum(abs(W(idx(j,:),j)-w_12)); */
			for (k = 0; k < p - 1; ++k) {
				if (sparse_W && a_(zero_ws, skip[k], j, p))
					continue;
				abs_change += fabs(a_(W, skip[k], j, p) - w_12[k]);
			}
			/* % update the column j and row j in W */
			/* W(idx(j,:),j) = w_12; */
			/* W(j,idx(j,:)) = w_12; */
			for (k = 0; k < p - 1; ++k) {
				if (sparse_W && a_(zero_ws, skip[k], j, p))
					continue;
				a_(W, skip[k], j, p) = w_12[k];
				a_(W, j, skip[k], p) = w_12[k];
			}
		}
		/* % calc mean change  */
		/* mean_abs_change = abs_change/(p^2-p); */
		mean_abs_change = abs_change * mean_abs_change_norm;
	}

	/* % calc precision matrix (see (7) and (8) in the glasso paper) */
	/* for j = 1:p */
	skip = skips; /* indexed by j */
	for (j = 0; j < p; ++j, skip += (p - 1)) {
		/* P(j,j) = 1 / (W(j,j) - W(j,idx(j,:))*beta_storage(j,:)'); */
		/* make oder of additions identical to Matlab */
		double P_jj = 0;
		for (k = 0; k < p - 1; ++k)
			P_jj += a_(W, j, skip[k], p) * b_(beta_storage, j, k, p);
		P_jj = 1 / (a_(W, j, j, p) - P_jj);
		a_(P, j, j, p) = P_jj;
		/* 	P(idx(j,:),j) = -beta_storage(j,:)*P(j,j); */
		/* 	P(j,idx(j,:)) = -beta_storage(j,:)*P(j,j); */
		for (k = 0; k < p - 1; ++k) {
			a_(P, skip[k], j, p) = -b_(beta_storage, j, k, p) * P_jj;
			a_(P, j, skip[k], p) = -b_(beta_storage, j, k, p) * P_jj;
		}
	}
	/* cleanup stuff */
	mxFree(skips);
	mxDestroyArray(w_12_array);
	mxDestroyArray(beta_storage_array);
	if (sparse_W)
		mxDestroyArray(zero_ws_array);
	if (sparse_P)
		mxDestroyArray(zero_betas_array);
	/* set return values */
	plhs[0] = W_array;
	plhs[1] = P_array;

	/* mexPrintf("cglasso bids farewell.\n"); */
}
