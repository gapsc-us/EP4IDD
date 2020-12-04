function [x_decod,prob_b]=PKSEPalg(A,complexFlag,sigma,r,h,pui_ldpc,nturbo,EPiterations)
    
% Function: PKSEPalg
%
% [x_decod,prob_b]=PKSEPalg(A,complexFlag,sigma,r,h,pui_ldpc,nturbo)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 25/07/2018
%
% Description: This function runs the P-KSEP [Santos18c] for turbo equalization. 
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
% r is the received signal
% h is the channel
% pui_ldpc is the probability of symbols at the output of the channel
% decoder
% nturbo is the number of the current turbo iteration
% EPiterations is parameter S, default value is 3
%
% Output: 
% x_decod is the estimation of the transmitted symbols
% prob_b is the (posterior or extrinsic) probability for each symbol
%
% References: 
% [Santos18c] = I. Santos, J. J. Murillo-Fuentes, and E. Arias-de-Reyna, 
% ?Equalization with expectation propagation at smoothing level, ? IEEE 
% Trans. On Sig. Proc., 2018. Under review. [Online]. Available: 
% https://arxiv.org/abs/1809.00806

flagextrinsic=1; % to 1 if the probability sent to the decoder is the extrinsic (cavity). to 0 if it is the posterior. 

M=length(A);
ntaps=length(h);
window=ntaps;

output_length=length(r);
input_length=output_length-ntaps+1;
AlphabetMat=(diag(A)*ones(M,input_length)).'; 

Hall=obth(h,window+1);
Hmat=Hall(ntaps:window,1:window);

varmin=1e-3;
A0=0;

%% EP parameters in [Santos18c]
if nargin==7
    L_EP=3;
else
    L_EP=EPiterations;
end
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
    sigma=sigma*sqrt(2);
    varLDPC=varLDPC*2;
end

nub=zeros(input_length,1);
varb=zeros(input_length,1);

%% EP in the outer loop (after the channel decoder)

for l=1:L_EP+1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Backwards %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    h_flip=fliplr(h);
    
    % Mean and variance of q in each iteration
    nu_qback=zeros(window,output_length);
    Cinv_qback=zeros(window,window,output_length);
    
    % Hacemos el backwards como un forward dándole la vuelta a los taps, u e y.
    hb=fliplr(h);
    hb_flip=fliplr(hb);
    rb=flipud(r);
    
    meanLDPC_b=flipud(meanLDPC);
    varLDPC_b=flipud(varLDPC);
    
    % Initialization
    meanWin=A0*ones(window,1);
    varWin=varmin*ones(window,1);
    C_q=varmin*eye(window);
    nu_q=A0*ones(window,1);
    gamma=meanWin./varWin;
    lambda=1./varWin;
    
    for kk=1:output_length
        if kk>input_length
            meanWin=[meanWin(2:end,:); A0];
            varWin=[varWin(2:end,:); varmin];
        else
            meanWin=[meanWin(2:end,:); meanLDPC_b(kk)];
            varWin=[varWin(2:end,:); varLDPC_b(kk)];
        end
        
        Cq_inv=inv(C_q(2:end,2:end));
        C_inv=sigma^(-2)*[zeros(window-ntaps,window);zeros(ntaps,window-ntaps) hb_flip'*hb_flip]+[Cq_inv-diag(lambda(2:end)) zeros(window-1,1); zeros(1,window)];
        zCinv=sigma^(-2)*[zeros(window-ntaps,1); hb_flip'*rb(kk)]+[Cq_inv*nu_q(2:end,1)-gamma(2:end,1);0];
        
        gamma=meanWin./varWin;
        lambda=1./varWin;
        
        C_q=inv(C_inv+diag(lambda));
        nu_q=C_q*(zCinv+gamma);
            
        nu_qback(:,output_length-kk+1)=nu_q(window:-1:1);
        Cinv_qback(:,:,output_length-kk+1)=inv(C_q(window:-1:1,window:-1:1));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Forwards %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    meanWin=A0*ones(window,1);
    varWin=varmin*ones(window,1);
    C_q=varmin*eye(window);
    nu_q=A0*ones(window,1);
    gamma=meanWin./varWin;
    lambda=1./varWin;
    
    for kk=1:input_length
        meanWin=[meanWin(2:end,:); meanLDPC(kk)];
        varWin=[varWin(2:end,:); varLDPC(kk)];
        
        Cq_inv=inv(C_q(2:end,2:end));
        C_inv=sigma^(-2)*[zeros(window-ntaps,window);zeros(ntaps,window-ntaps) h_flip'*h_flip]+[Cq_inv-diag(lambda(2:end)) zeros(window-1,1); zeros(1,window)];
        zCinv=sigma^(-2)*[zeros(window-ntaps,1); h_flip'*r(kk)]+[Cq_inv*nu_q(2:end,1)-gamma(2:end,1);0];
        
        gamma=meanWin./varWin;
        lambda=1./varWin;
        
        C_q=inv(C_inv+diag(lambda));
        nu_q=C_q*(zCinv+gamma);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Smoothing %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Denominator
        Cinv_p=sigma^(-2)*Hmat'*Hmat;
        zCinv_p=sigma^(-2)*Hmat'*r(kk);
        Covinv_den=Cinv_p+diag(lambda);
        Covinvnu_den=zCinv_p+gamma;
        
        % Numerator
        C_qinvf=inv(C_q);
        Covinv_num=C_qinvf+Cinv_qback(:,:,kk);
        Covinvnu_num=C_qinvf*nu_q+Cinv_qback(:,:,kk)*nu_qback(:,kk);
        
        % Smoothing
        Cinv_qall=Covinv_num-Covinv_den;
        Cinvnu_qall=Covinvnu_num-Covinvnu_den;
        
        C_qEP=inv(Cinv_qall);
        var_qi=diag(C_qEP);
        nu_qi=C_qEP*(Cinvnu_qall);
        
        % Extrinsic
        Ve=var_qi./(1-var_qi.*lambda);
        me=Ve.*(nu_qi./var_qi-gamma);
        
        if flagextrinsic
            nub(kk)=me(window);
            varb(kk)=Ve(window);
        else
            nub(kk)=nu_qi(window);
            varb(kk)=var_qi(window);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% EP %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        if l<=L_EP
    
            % Compute cavity marginal
            if complexFlag
                cavity_i=normpdfcomplex(A,me(window),sqrt(Ve(window)));
            else
                cavity_i=normpdf(A,me(window),sqrt(Ve(window)));
            end
            cavity_i=max(cavity_i,minprob); % to avoid NaN when normalizing
            
            % Normalize distribution
            K=sum(cavity_i);
            cavity_i=cavity_i./K;
            
            posterior=cavity_i.*pui_ldpc(kk,:);
            posterior=max(posterior,minprob); % to avoid NaN when normalizing
            K=sum(posterior);
            posterior=posterior./K;
            
            % Mean (nu_pi) and variance (var_pi) of distribution (33)
            nu_pi=sum((diag(A)*transpose(posterior)).',2);
            var_pi=sum(conj(A-nu_pi).*(A-nu_pi).*posterior,2);
            var_pi=max(var_pi,epsilon);
            
            % Compute new values for gamma and GAMMA
            GAMMA_aux=(1./var_pi)-(1./Ve(window));
            gamma_aux=(nu_pi./var_pi)-(me(window)./Ve(window));
            
            GAMMAupdated=beta*GAMMA_aux+(1-beta)*lambda(window);
            gammaupdated=beta*gamma_aux+(1-beta)*gamma(window);
            
            if GAMMAupdated>0 && sum(isfinite(cavity_i))==M
                meanLDPC(kk)=gammaupdated/GAMMAupdated;
                varLDPC(kk)=1/GAMMAupdated;
            end
        end
    end
end

%% Probabilities at the output of the equalizer

if complexFlag
    prob_b=normpdfcomplex(AlphabetMat,nub*ones(1,M),sqrt(varb)*ones(1,M)).';
else
    prob_b=normpdf(AlphabetMat,nub*ones(1,M),sqrt(varb)*ones(1,M)).';
end
prob_b=prob_b./transpose(diag(sum(prob_b))*ones(input_length,M));


[~,ind]=max(prob_b);
x_decod=A(ind).';
    