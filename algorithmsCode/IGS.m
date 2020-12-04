function [invW,x]=IGS(W,b,K)

% IGS  Improved Gauss-Seidel (IGS) detection for MIMO
%   invW,invWb=IGS(W,b) obtains an approximation to the inv(W) based on the
%   Neumman serie of first order and also obtains an approximation for
%   x=inv(W)*b. The n-by-n matrix W must be symmetric and 
%   the column vector b must have length n.

% INPUTS
% W  : matrix n x n
% b  : vector to solve against n x 1
% K  : number of iterations of the algorithm
%
% OUTPUT
% invW : approximation to inv(W)
% x : approximation to x=inv(W)*b
%
% Copyright (c) Irene Santos, 2019.
%
% Reference: C. Zhang, Z. Wu, C. Studer, Z. Zhang, X. You, "Efficient 
% Soft-Output Gauss-Seidel Data Detector for Massive MIMO Systems", 2018. 
% https://arxiv.org/abs/1804.06737

if nargin<3 || isempty(K)
    K=3;
end

D=diag(diag(W));
invD=diag(1./diag(D));
E=W-D;
L=tril(E);

invW=invD-invD*E*invD; % approximation to the inverse based on the neumann serie of first order
% Initialization of x
x=invW*b;

for k=1:K
    x=(D+L)\(b-L'*x);
end