function [x_decod,prob_b]=BEPalg(A,complexFlag,sigma,y,H,pui_ldpc,S_EP)
%
% Function: BEPalg
%
% [x_decod,prob_b]=BEPalg(A,complexFlag,sigma,y,H,pui_ldpc, S_EP)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 20/02/2018
%
% Description: This function runs the BEP [Santos16] for turbo equalization. 
% It returns the estimated symbols (x_decod) and the probability for each 
% symbol (prob_b). If flagextrinsic=1, then prob_b is the extrinsic 
% distribution. If flagextrinsic=0, then prob_b is the posterior
% distribution. 
% 
% Inputs: 
%  A is the set of symbols
%   complexFlag indicates if the symbols are complex or real (1-complex,
%   0-real)
%  sigma is the standard deviation of the noise
%  y is the received signal
%  H is the channel matrix
%  pui_ldpc is the probability of symbols at the output of the channel
%   decoder
%  S_EP: number of iterations of the (inner) EP, default value is 10
%
% Output: 
%  x_decod is the estimation of the transmitted symbols
%  prob_b is the (posterior or extrinsic) probability for each symbol
%
% References: 
% [Santos16] = I. Santos, J. J. Murillo-Fuentes, R. Boloix-Tortosa, E. 
% Arias-de-Reyna, and P. M. Olmos, ?Expectation propagation as turbo 
% equalizer in ISI channels,? IEEE Transactions on Communications, vol. 65, 
% no. 1, pp. 360-370, Jan 2017. 
% [Cespedes18] = J. Céspedes, P. M. Olmos, M. Sánchez-Fernández and F. 
% Pérez-Cruz, ?Probabilistic MIMO symbol detection with expectation 
% consistency approximate inference,? in IEEE Trans. On Vehicular Tech., 
% vol. 67, no. 4, pp. 3481-3494, April 2018.

input_length=size(H,2);
M=length(A);

flagextrinsic=1; % to 1 if the probability sent to the decoder is the extrinsic (cavity). to 0 if it is the posterior. 

AlphabetMat=(diag(A)*ones(M,input_length)).'; 

%% EP parameters in [Santos16]
if nargin<7
S_EP=10; % in [Santos16]
end

beta=0.1;
if isempty(pui_ldpc) % standalone equalization
    if M>=64
        eps0=0.9;
    else
        eps0=-1; % it means that the value of eps is computed along the l-th iterations of the EP algorithm
    end
else % turbo equalization
    eps0=1e-8;
end
epsilon=eps0;

% For running the EPD in [Cespedes18], we just have to set the EP
% parameters to: 
% beta=0.95;
% eps0=-1;
% L_EP=10;

minprob=1e-8; % minimum value for the probabilities

%% Initialization of priors
if isempty(pui_ldpc)
    pui_ldpc=1/M*ones(input_length,M);
end

[meanLDPC,varLDPC]=meanvar(pui_ldpc.',A);
if complexFlag
    % The complex system can be simulated by: 
    % - Considering the real and imaginary parts separately 
    % - Formulating the algorithms in complex formulation
    % The following is an adjustment to let the results of both
    % possibilities match: 
    sigma=sigma*sqrt(2);
    varLDPC=varLDPC*2;
end

GAMMA=1./varLDPC;
gamma=meanLDPC./varLDPC;

%% Initial iteration --> LMMSE

% meand and variance vectors of all marginal qi(x) at iteration l=0
C_q=inv(sigma^(-2)*H'*H+diag(GAMMA));
var_qi=diag(C_q);
nu_qi=C_q*(sigma^(-2)*H'*y+gamma);

% Mean (ti) and variance (hi_2) of each cavity marginal
hi_2=var_qi./(1-var_qi.*GAMMA);
ti=hi_2.*(nu_qi./var_qi-gamma);

%% EP in the inner loop

for l=1:S_EP
    
    if eps0<0
        epsilon=max([2.^(-max([l-5, 1])),1e-8]);
    end
      
    % Compute cavity marginal for each xi
    if complexFlag
        cavity_i=normpdfcomplex(AlphabetMat,ti*ones(1,M),sqrt(hi_2)*ones(1,M));
    else
        cavity_i=normpdf(AlphabetMat,ti*ones(1,M),sqrt(hi_2)*ones(1,M));
    end
    
    cavity_i=max(cavity_i,minprob); % to avoid NaN when normalizing
    
    % Normalize distribution (33)
    K=sum(cavity_i,2);
    cavity_i=cavity_i./(diag(K)*ones(input_length,M));
    
    % Mean (nu_pi) and variance (var_pi) of distribution (33)
    nu_pi=sum((diag(A)*transpose(cavity_i)).',2);
    var_pi=sum(conj(AlphabetMat-diag(nu_pi)*ones(input_length,M)).*(AlphabetMat-diag(nu_pi)*ones(input_length,M)).*cavity_i,2);
    var_pi=max([var_pi,epsilon*ones(input_length,1)],[],2);
    
    % Compute new values for gamma and GAMMA
    GAMMA_aux=(1./var_pi)-(1./hi_2);
    gamma_aux=(nu_pi./var_pi)-(ti./hi_2);
    
    % damping
    GAMMAupdated=beta*GAMMA_aux+(1-beta)*GAMMA;
    gammaupdated=beta*gamma_aux+(1-beta)*gamma;
    
    % negative variances?
    indexNeg=find(GAMMAupdated>0 & sum(isfinite(cavity_i),2)==M);
    GAMMA(indexNeg)=GAMMAupdated(indexNeg);
    gamma(indexNeg)=gammaupdated(indexNeg);
    
    % meand and variance vectors of all marginal qi(x) at next
    % iteration l
    C_q=inv(sigma^(-2)*H'*H+diag(GAMMA));
    var_qi=diag(C_q);
    nu_qi=C_q*(sigma^(-2)*H'*y+gamma);
    
    % Mean (ti) and variance (hi_2) of each cavity marginal
    hi_2=var_qi./(1-var_qi.*GAMMA);
    ti=hi_2.*(nu_qi./var_qi-gamma);
end

%% Probabilities at the output of the equalizer
% Decoder
if flagextrinsic
    if complexFlag
        prob_b=normpdfcomplex(AlphabetMat,ti*ones(1,M),sqrt(hi_2)*ones(1,M)).';
    else
        prob_b=normpdf(AlphabetMat,ti*ones(1,M),sqrt(hi_2)*ones(1,M)).';
    end
else
    if complexFlag
        prob_b=normpdfcomplex(AlphabetMat,nu_qi*ones(1,M),sqrt(var_qi)*ones(1,M)).';
    else
        prob_b=normpdf(AlphabetMat,nu_qi*ones(1,M),sqrt(var_qi)*ones(1,M)).';
    end
end

prob_b=prob_b./transpose(diag(sum(prob_b))*ones(input_length,M));
[~,ind]=max(prob_b);
x_decod=A(ind).';
