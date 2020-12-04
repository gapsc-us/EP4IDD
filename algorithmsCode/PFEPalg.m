function [x_decod,prob_b,LLRs,hi_2,ti]=PFEPalg(A,complexFlag,sigma,y,h,H,pui_ldpc,nturbo)

% Function: PFEPalg
%
% [x_decod,prob_b,LLRs,hi_2,ti]=PFEPalg(A,complexFlag,sigma,y,h,H,pui_ldpc,nturbo)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 27/02/2017
%
% Description: This function runs the P-FEP [Santos18] for turbo equalization. 
% It returns the estimated symbols (x_decod) and the probability for each 
% symbol (prob_b). If flagextrinsic=1, then prob_b is the extrinsic 
% distribution. If flagextrinsic=0, then prob_b is the posterior
% distribution. 
% 
% Inputs: 
% A is the set of symbols
% complexFlag indicates if the symbols are complex or real (1-complex,
% 0-real)
% sigma is the standard deviation of the noise
% y is the received signal
% h is the channel
% H is the channel matrix
% pui_ldpc is the probability of symbols at the output of the channel
% decoder
% nturbo is the number of the current turbo iteration
%
% Output: 
% x_decod is the estimation of the transmitted symbols
% prob_b is the (posterior or extrinsic) probability for each symbol
% LLRs is the loglikehood ratio at the output of the equalizer
% ti is the mean of the extrinsic distribution at the output of the equalizer
% hi_2 is the variance of the extrinsic distribution at the output of the 
% equalizer
%
% References: 
% [Santos18] = I. Santos, J. J. Murillo-Fuentes, E. Arias-de-Reyna, and P. 
% M. Olmos, ?Turbo EP-based equalization: A filter-type implementation,? 
% IEEE Trans. On Comm., vol. 66, no. 9, pp. 4259-4270, Sept. 2018

input_length=size(H,2);
M=length(A);
energy=sum(abs(A).^2)/M;

% Limits of the window in Tuchler
ntaps=length(h);
N1=2*ntaps;N2=ntaps+1; 

flagextrinsic=1; % to 1 if the probability sent to the decoder is the extrinsic (cavity). to 0 if it is the posterior. 

AlphabetMat=(diag(A)*ones(M,input_length)).'; 

%% EP parameters in [Santos18]
L_EP=3;
beta=min(exp(nturbo/1.5)/10,0.7);
epsilon=1e-8;

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
    energy=energy*2;
    sigma=sigma*sqrt(2);
    varLDPC=varLDPC*2;
end

GAMMA=1./varLDPC;
gamma=meanLDPC./varLDPC;

%% Initial iteration --> LMMSE (Wiener)

[ti,LLRs,hi_2]=mmseTuchler(y(1:input_length),N1,N2,h,sigma,meanLDPC,varLDPC,energy);

%% EP in the inner loop

for l=1:L_EP
    
    % Compute cavity marginal for each xi
    if complexFlag
        cavity_i=real(normpdfcomplex(AlphabetMat,ti*ones(1,M),sqrt(hi_2)*ones(1,M)));
    else
        cavity_i=normpdf(AlphabetMat,ti*ones(1,M),sqrt(hi_2)*ones(1,M));
    end
    cavity_i=max(cavity_i,minprob); % to avoid NaN when normalizing
    
    % Normalize distribution (33)
    K=sum(cavity_i,2);
    cavity_i=cavity_i./(diag(K)*ones(input_length,M));
    
    posterior=cavity_i.*pui_ldpc;
    posterior=max(posterior,minprob); % to avoid NaN when normalizing
    K=sum(posterior,2);
    posterior=posterior./(diag(K)*ones(input_length,M));
    
    % Mean (nu_pi) and variance (var_pi) of distribution (33)
    nu_pi=sum((diag(A)*posterior').',2);
    var_pi=sum(conj(AlphabetMat-diag(nu_pi)*ones(input_length,M)).*(AlphabetMat-diag(nu_pi)*ones(input_length,M)).*posterior,2);
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
    
    %con Tuchler    
    varLDPC=1./GAMMA;       
    meanLDPC=gamma./GAMMA;   
    [ti,LLRs,hi_2]=mmseTuchler(y(1:input_length),N1,N2,h,sigma,meanLDPC,varLDPC,energy);
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

prob_b=real(prob_b./(diag(sum(prob_b))*ones(input_length,M)).');
[~,ind]=max(prob_b);
x_decod=A(ind).';