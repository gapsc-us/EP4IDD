function [x_decod,prob_b,me,Ve]=MMSEalgAproxNeuman(A,complexFlag,sigma,y,H,pui_ldpc)

% Function: MMSEalgAproxNeuman
%
% [x_decod,prob_b,me,Ve]=MMSEalgAproxNeuman(A,complexFlag,sigma,y,H,pui_ldpc)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 7/03/2017
%
% Description: This function runs an approximation to the LMMSE [Muranov10] 
% for turbo equalization based on the Neumann Series expansion [Fang16]. 
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
% H is the channel matrix
% pui_ldpc is the probability of symbols at the output of the channel
% decoder
%
% Output: 
% x_decod is the estimation of the transmitted symbols
% prob_b is the (posterior or extrinsic) probability for each symbol
% me is the mean of the extrinsic distribution at the output of the equalizer
% Ve is the variance of the extrinsic distribution at the output of the 
% equalizer
%
% References: 
% [Muranov10] = K. Muranov, "Survey of MMSE channel equalizers," Tech. 
% Rep., University of Illinois, Chicago. 
%
% [Fang16] = L. Fang, L. Xu and D. D. Huang, "Low Complexity Iterative 
% MMSE-PIC Detection for Medium-Size Massive MIMO," in IEEE Wireless 
% Communications Letters, vol. 5, no. 1, pp. 108-111, Feb. 2016.

input_length=size(H,2);
M=length(A);
AlphabetMat=(diag(A)*ones(M,input_length)).'; 

flagextrinsic=1; % to 1 if the probability sent to the decoder is the extrinsic (cavity). to 0 if it is the posterior. 

%% Parameters of RCNSE
Knse=3;

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

%% LMMSE (Block)

% meand and variance vectors of all marginal qi(x) at iteration l=0
[z,varC]=RCNeumannSerie(meanLDPC,varLDPC,H,y,sigma^2,Knse);

%% Probabilities at the output of the equalizer

% Decoder

if flagextrinsic
    Ve=1./(1./varC-GAMMA);
    me=Ve.*(z./varC-gamma);
    
    if complexFlag
        prob_b=normpdfcomplex(AlphabetMat,me*ones(1,M),sqrt(Ve)*ones(1,M))';
    else
        prob_b=normpdf(AlphabetMat,me*ones(1,M),sqrt(Ve)*ones(1,M))';
    end
else
    if complexFlag
        prob_b=normpdfcomplex(AlphabetMat,z*ones(1,M),sqrt(diag(C))*ones(1,M))';
    else
        prob_b=normpdf(AlphabetMat,z*ones(1,M),sqrt(diag(C))*ones(1,M))';
    end
end

prob_b=prob_b./transpose(diag(sum(prob_b))*ones(input_length,M));
[~,ind]=max(prob_b);
x_decod=A(ind).';

        
