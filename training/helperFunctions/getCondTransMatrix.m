function matrix = getCondTransMatrix(mu,Prec,numRows,doTriu);

% getCondTransMatrix(mu,Prec,numRows,options) - evaluates the distribution p(a|b) at discrete positions 1:numRows 
%
% Syntax:
%   matrix = getCondTransMatrix(mu,Prec,numRows,options)
%
% Inputs:
%   mu      - [array](1x2) the mean of a and b
%   Prec    - [matrix](2x2) precision matrix for two neighboring boundaries
%   numRows - [int] number of rows within the B-Scan
%
% Outputs:
%   matrix - [matrix](numRowsxnumRows) transition matrix for p(a|b) 
%
% See also: trainShape, setShapeDefaults

% Author: Fabian Rathke
% email: frathke@googlemail.com
% Website: https://github.com/FabianRathke/octSegmentation
% Last Revision: 06-Dec-2013

positions = 1:numRows;
% mean of a conditioned on b
mu_a_b = mu(2) - Prec(2,2)^-1*Prec(2,1)*(positions-mu(1));

A = positions;
factor = 1/(2*pi*inv(Prec(2,2)))^0.5;
n = length(positions);
tmp = reshape(bsxfun(@minus,A',mu_a_b),1,n^2);
if nargin > 3 & ~doTriu
	matrix = factor*reshape(exp(-0.5*Prec(2,2)*(tmp.*tmp)),n,n)';
else
	matrix = triu(factor*reshape(exp(-0.5*Prec(2,2)*(tmp.*tmp)),n,n)',0);
end

% introduce sparsity
matrix(matrix<(10^-20)) = 0;
