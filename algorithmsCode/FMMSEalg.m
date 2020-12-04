function [x_decod,prob_b,LLRs,hi_2,ti]=FMMSEalg(A,complexFlag,sigma,y,h,H,pui_ldpc)

% Function: FMMSEalg
%
% [x_decod,prob_b,LLRs,hi_2,ti]=FMMSEalg(A,complexFlag,sigma,y,h,H,pui_ldpc)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 30/12/2016
%
% Description: This function runs the MMSE Wiener [Tuchler02a,Tuchler02b,
% Tuchler11] for turbo equalization. 
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
% [Tuchler02a] = M. Tüchler, R. Koetter, and A. Singer, ?Turbo 
% equalization: principles and new results,? IEEE Trans. On Communications, 
% vol. 50, no. 5, pp. 754-767, May 2002.
% [Tuchler02b] = M. Tüchler, A. Singer, and R. Koetter, ?Minimum mean 
% squared error equalization using a priori information, ? IEEE Trans. On 
% Signal Processing, vol. 50, no 3, pp. 673-683, Mar 2002.
% [Tuchler11] = M. Tüchler, and A. Singer, ?Turbo equalization: An 
% overview,? IEEE Trans. On Information Theory, vol. 57, no. 2, pp. 920-952, 
% Feb 2011. 

input_length=size(H,2);
M=length(A);
energy=sum(abs(A).^2)/M;

% Limits of the window in Tuchler
ntaps=length(h);
N1=2*ntaps;N2=ntaps+1; 

flagextrinsic=1; % to 1 if the probability sent to the decoder is the extrinsic (cavity). to 0 if it is the posterior. 

AlphabetMat=(diag(A)*ones(M,input_length)).'; 

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

%% LMMSE (Wiener)

[ti,LLRs,hi_2]=mmseTuchler(y(1:input_length),N1,N2,h,sigma,meanLDPC,varLDPC,energy);

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