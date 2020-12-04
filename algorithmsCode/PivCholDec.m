function Lk=PivCholDec(A,k)

% PivCholDec  Pivoted Cholesky factorizacion 
%   Lk=PivCholDec(A,k) obtains a low-rank approximation of a positive
%   definite matrix A by means of the pivoted cholesky factorization
%
% INPUTS
% A  : matrix n x n
% k  : rank of the approximation
%
% OUTPUT
% Lk : low-rank approximation, i.e., A \approx Lk*Lk'
%
% Copyright (c) Irene Santos, 2019.
%
% References: 
% [Harbrecht12] = Harbrecht, H., Peters, M., Schneider, R. (2012). On the 
% low-rank approximation by the pivoted Cholesky decomposition. Applied 
% numerical mathematics, 62(4), 428-440.
%
% [Gardner18] = J. R. Gardner, G. Pleiss, D. Bindel, K. Q. Weinberger, A.
% G. Wilson, "GPyTorch: Blackbox Matrix-Matrix Gaussian Process Inference
% with GPU Acceleration" 2018, [Online] https://arxiv.org/abs/1809.11165

n=size(A,1);

if nargin<2 || isempty(k)
    k=n;
end
k=min(k,rank(A));

P=zeros(n,n,k);
q=zeros(n,1,k);
MulP=eye(n);
Lk=[];

Sk=A;

for ii=1:k
    % maximum diagonal element of Sk
    diagA=diag(Sk);
    [~,index]=max(diagA);
    
    % Permuting matrix
    Phat=eye(size(Sk));
    aux=Phat(:,1);
    Phat(:,1)=Phat(:,index);
    Phat(:,index)=aux;
    P(:,:,ii)=eye(n);
    P(ii:n,ii:n,ii)=Phat;
    
    % Permuting Sk
    Sk=P(ii:n,ii:n,ii)*Sk*P(ii:n,ii:n,ii);
    
    % Cholesky decomposition
    a11=Sk(1,1);
    b=Sk(2:end,1);
    A22=Sk(2:end,2:end);
    q(ii:n,:,ii)=[a11;b]/sqrt(a11);
    Sk=A22-b*b'/a11;
    
    % multiplication of permutation matrices
    MulP=MulP*P(:,:,ii);
    Lk=[Lk MulP*q(:,:,ii)];      
end

