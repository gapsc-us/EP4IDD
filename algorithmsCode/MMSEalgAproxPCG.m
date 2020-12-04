function [x_decod,prob_b,me,Ve]=MMSEalgAproxPCG(A,complexFlag,sigma,y,H,pui_ldpc)

% Function: MMSEalgAproxPCG
%
% [x_decod,prob_b,me,Ve]=MMSEalgAproxPCG(A,complexFlag,sigma,y,H,pui_ldpc)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 7/03/2017
%
% Description: This function runs an approximation to the LMMSE [Muranov10] 
% for turbo equalization based on the conjugate gradient (CG) [Yin14]. A
% preconditioning matrix can be also included, yielding the preconditioned 
% conjugate gradient (PCG) [Gardner18]. 
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
% [Yin14] = B. Yin, M. Wu, J. R. Cavallaro, and C. Studer, "Conjugate 
% gradient- based soft-output detection and precoding in massive MIMO 
% systems," in IEEE Global Comm. Conf., Dec 2014, pp. 3696-3701.
%
% [Gardner18] = J. R. Gardner, G. Pleiss, D. Bindel, K. Q. Weinberger, A.
% G. Wilson, "GPyTorch: Blackbox Matrix-Matrix Gaussian Process Inference
% with GPU Acceleration" 2018, [Online] https://arxiv.org/abs/1809.11165

input_length=size(H,2);
output_length=size(H,1);
M=length(A);
AlphabetMat=(diag(A)*ones(M,input_length)).'; 

flagextrinsic=1; % to 1 if the probability sent to the decoder is the extrinsic (cavity). to 0 if it is the posterior. 

%% Parameters of CG
Kcg=4; % number of iterations of the CG
flagprecond=0; % To 1 if we use a preconditioner matrix. 0 in other case

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

% Matrix to invert: Nt x Nt
G=sigma^(-2)*(H'*H);
K=G+diag(GAMMA);
b=H'*y-H'*(H*(gamma./GAMMA));
%b=sigma^(-2)*H'*y+gamma;
[invK,~]=IGS(K,b,0); % Neuman serie of order 1 for the variance
varC=diag(invK);
if flagprecond
    Lk=PivCholDec(G,10); % If the second input parameter is size(G,1), it is equivalent to Lk=chol(G, 'lower');
    Pk=Lk*Lk';
    Phat = Pk+diag(GAMMA); % preconditioning
else
    Phat=eye(input_length); % no preconditioner
end
[z,~,~,iter] = pcg(K,b,1e-6,Kcg,Phat,Phat');
z=sigma^(-2)*z+gamma./GAMMA;

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

        
