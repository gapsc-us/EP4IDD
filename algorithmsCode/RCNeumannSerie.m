function [m_posterior,v_posterior]=RCNeumannSerie(m_priori,v_priori,H,y,sigma2,L)

% RCNeumannSerie  Reduced complexity Neuman Series expansion for a
% turbo detector with received signal y, matrix channel H, noise variance
% sigma2 and a priori mean (m_priori) and variance (v_priori) on the
% transmitted symbols. 
%   [m_posterior,v_posterior]=RCNeumannSerie(m,V,H,sigma2,L) computes the a 
%   posteriori mean and variance for a MMSE detector with a priori 
%   information about the transmitted symbols (i.e., turbo MMSE). It
%   obtains an approximation to: 
%   V_posterior = inv(diag(1./v_priori)+H'*H/sigma2);
%   m_posterior = m_priori+V_posterior/sigma2*(H'*y-H'*H*m_priori)
%
% INPUTS
% m_priori  : vector Nt x 1 with the a priori means
% v_priori  : vector Nt x 1 with the a priori variances
% H  : matrix Nr x Nt with the channel weights
% y  : vector Nr x 1 with the received signal
% sigma2 : noise variance
% L  : number of iterations of the algorithm
%
% OUTPUT
% m_posterior : vector Nt x 1 with the a priori means
% v_posterior : vector Nt x 1 with the a priori variances
%
% Copyright (c) Irene Santos, 2019. 
%
% Reference: L. Fang, L. Xu and D. D. Huang, "Low Complexity Iterative 
% MMSE-PIC Detection for Medium-Size Massive MIMO," in IEEE Wireless 
% Communications Letters, vol. 5, no. 1, pp. 108-111, Feb. 2016.

if isempty(L) || ~exist('L','var')
    L=3;
end

invV=diag(1./v_priori);

% Exact mean and variance
%V_posterior = inv(invV+H'*H/sigma2);
% v_posterior=diag(V_posterior);
%m_posterior = m_priori+V_posterior/sigma2*(H'*y-H'*H*m_priori);

% Calculate a posteriori mean v_posterior
D=invV+diag(diag(H'*H)/sigma2); % We only need the diagonal elements of H'*H
invD=diag(1./diag(D)); % We use the diagonal nature of D to invert efficiently
v=invD*(H'*y-H'*(H*m_priori)); % We do not have to perform a matrix-matrix multiplication if we first compute the vector H*m and then multiply this vector by H'
s=v;
for ii=1:L
    v=v-invD*(invV*v+H'*(H*v/sigma2)); % Again, we first compute the vector H*v and then multiply it by H' rather than performing the matrix-matrix multiplication H'*H
    s=s+v;
end
m_posterior=m_priori+s/sigma2;

% Approximate the diagonal elements of V_posterior
v_posterior=diag(invD);

