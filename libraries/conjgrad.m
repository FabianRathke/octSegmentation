function x = conjgrad(A,b,x)
% conjgrad - solves the system of linear equations Ax = b via conjugate gradient (code example from http://en.wikipedia.org/wiki/Conjugate_gradient_method)
% 
% Syntax:
%   x = conjgrad(A,b,x)
%
% Inputs:
%   A - [matrix]
%   b - [array] 
%   x - [array] initial solution; can be a vector of all zeros
%
% Outputs:
%   x - [array] the final solution vector 
%
% See also: trainAppearance

% Author: Fabian Rathke (Not my code: http://en.wikipedia.org/wiki/Conjugate_gradient_method)
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 12-Dec-2013

r=b-A*x;
p=r;
rsold=r'*r;

for i=1:10000000
	Ap=A*p;
	alpha=rsold/(p'*Ap);
	x=x+alpha*p;
	r=r-alpha*Ap;
	rsnew=r'*r;
	if sqrt(rsnew)<1e-10
		  break;
	end
	p=r+rsnew/rsold*p;
	rsold=rsnew;
end
