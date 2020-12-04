function y = mvnpdfcomplex(X, MU, SIGMA)

% Function: mvnpdfcomplex
%
% y = mvnpdfcomplex(X, MU, SIGMA)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 7/03/2017
%
% Description: This function computes the multivariate normal probability 
% density function (pdf) with mean MU and covariance SIGMA, evaluated at 
% each row of X (where X is a matrix with complex numbers)
% 
% Inputs: 
% X is a N-by-D complex matrix where the rows correspond to points and
% columns correspond to variables. The pdf will be evaluated at these
% points. 
% MU is a N-by-D matrix with the complex mean
% SIGMA is a D-by-D matrix with the variances. 
%
% Output: 
% y is a N-by-1 vector with the pdf

% Get size of data.  Column vectors provisionally interpreted as multiple scalar data.
[N,D] = size(X);

% Checking if SIGMA is diagonal
if isdiag(SIGMA)
    sigmaIsDiag = true;
else
    sigmaIsDiag = false;
end

if sigmaIsDiag
    sigmaInv=diag(1./diag(SIGMA));
    detSigma=prod(diag(SIGMA));
else
    sigmaInv=SIGMA\eye(length(SIGMA));
    detSigma=det(SIGMA);
end

y=exp(-ctranspose(X-MU)*sigmaInv*(X-MU))./(pi^N*detSigma);

