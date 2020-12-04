function [x_decod,prob_b]=PKSEPalg_OptComp(A,complexFlag,sigma,r,h,pui_ldpc,nturbo,EPiterations)
    
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
% symbol (prob_b). Its computational complexity is optimal, i.e., O(L^2)
% since in the forward/bacward steps we do not invert matrices (it is
% solved with matrix-vector multiplications). In the smoothing step, the
% inverses only have to be computed each L symbols. 
% If flagextrinsic=1, then prob_b is the extrinsic distribution. If 
% flagextrinsic=0, then prob_b is the posterior distribution. 
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
% EPiterations is the parameter S, default value is 3
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
AlphabetMatW=(diag(A)*ones(M,window)).'; % Matrix WxL

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
    Cq_inv=1./varmin*eye(window);
    nu_q=A0*ones(window,1);
    
    for kk=1:output_length
        if kk>input_length
            meanWin=[meanWin(2:end,:); A0];
            varWin=[varWin(2:end,:); varmin];
        else
            meanWin=[meanWin(2:end,:); meanLDPC_b(kk)];
            varWin=[varWin(2:end,:); varLDPC_b(kk)];
        end
        
        gamma=meanWin./varWin;
        lambda=1./varWin;
        
        Cq_inv_marg=Cq_inv(2:end,2:end)-Cq_inv(2:end,1)*Cq_inv(1,2:end)/Cq_inv(1,1);
        Cq_inv=sigma^(-2)*(hb_flip'*hb_flip)+[Cq_inv_marg zeros(window-1,1); zeros(1,window-1) lambda(window)]; % Inverse of eq. (8)
        zCinv=sigma^(-2)*hb_flip'*rb(kk)+[Cq_inv_marg*nu_q(2:end,1);gamma(window)];     
        
        B=[C_q(2:end,2:end) zeros(window-1,1); zeros(1,window-1) varWin(window)];
        C_q=B-B*(hb_flip'*hb_flip)*B/(sigma^(2)+hb_flip*B*hb_flip'); % eq. (8)
        nu_q=C_q*zCinv; % eq. (7)
        
        nu_qback(:,output_length-kk+1)=nu_q(window:-1:1);
        Cinv_qback(:,:,output_length-kk+1)=Cq_inv(window:-1:1,window:-1:1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Forwards %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    meanWin=A0*ones(window,1);
    varWin=varmin*ones(window,1);
    C_q=varmin*eye(window);
    Cq_inv=1./varmin*eye(window);
    nu_q=A0*ones(window,1);
    
    for kk=1:input_length
        meanWin=[meanWin(2:end,:); meanLDPC(kk)];
        varWin=[varWin(2:end,:); varLDPC(kk)];
        
        gamma=meanWin./varWin;
        lambda=1./varWin;
        
        Cq_inv_marg=Cq_inv(2:end,2:end)-Cq_inv(2:end,1)*Cq_inv(1,2:end)/Cq_inv(1,1);
        Cq_inv=sigma^(-2)*(h_flip'*h_flip)+[Cq_inv_marg zeros(window-1,1); zeros(1,window-1) lambda(window)]; % Inverse of eq. (8)
        zCinv=sigma^(-2)*h_flip'*r(kk)+[Cq_inv_marg*nu_q(2:end,1);gamma(window)];
        
        B=[C_q(2:end,2:end) zeros(window-1,1); zeros(1,window-1) varWin(window)];
        C_q=B-B*(h_flip'*h_flip)*B/(sigma^(2)+h_flip*B*h_flip'); % eq. (8)
        nu_q=C_q*zCinv; % eq. (7)
              
        if mod(kk,window)==0 || kk==input_length
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%% Smoothing %%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Denominator
            Cinv_p=sigma^(-2)*(Hmat'*Hmat);
            zCinv_p=sigma^(-2)*Hmat'*r(kk);
            Covinv_den=Cinv_p+diag(lambda);
            Covinvnu_den=zCinv_p+gamma;
            
            % Numerator
            Covinv_num=Cq_inv+Cinv_qback(:,:,kk);
            Covinvnu_num=Cq_inv*nu_q+Cinv_qback(:,:,kk)*nu_qback(:,kk);
            
            % Smoothing
            Cinv_qall=Covinv_num-Covinv_den;
            Cinvnu_qall=Covinvnu_num-Covinvnu_den;
            
            C_qEP=inv(Cinv_qall); % eq. (12)
            var_qi=diag(C_qEP);
            nu_qi=C_qEP*(Cinvnu_qall); % eq. (11)
            
            % Extrinsic
            Ve=var_qi./(1-var_qi.*lambda);
            me=Ve.*(nu_qi./var_qi-gamma);
            
            if flagextrinsic
                nub(kk-ntaps+1:kk)=me;
                varb(kk-ntaps+1:kk)=Ve;
            else
                nub(kk-ntaps+1:kk)=nu_qi;
                varb(kk-ntaps+1:kk)=var_qi;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%% EP %%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            if l<=L_EP
                % Compute cavity marginal
                if complexFlag
                    cavity_i=normpdfcomplex(AlphabetMatW,me*ones(1,M),sqrt(Ve)*ones(1,M));
                else
                    cavity_i=normpdf(AlphabetMatW,me*ones(1,M),sqrt(Ve)*ones(1,M));
                end
                cavity_i=max(cavity_i,minprob); % to avoid NaN when normalizing
                
                % Normalize distribution
                K=sum(cavity_i,2);
                cavity_i=cavity_i./(diag(K)*ones(window,M));
                
                % discrete posterior                
                posterior=cavity_i.*pui_ldpc(kk-ntaps+1:kk,:);
                posterior=max(posterior,minprob); % to avoid NaN when normalizing
                K=sum(posterior,2);
                posterior=posterior./(diag(K)*ones(window,M));
                
                % Mean (nu_pi) and variance (var_pi) of distribution (33)
                nu_pi=sum((diag(A)*transpose(posterior)).',2);
                var_pi=sum(conj(AlphabetMatW-diag(nu_pi)*ones(window,M)).*(AlphabetMatW-diag(nu_pi)*ones(window,M)).*posterior,2);
                var_pi=max([var_pi,epsilon*ones(window,1)],[],2);
                
                % Compute new values for gamma and GAMMA
                GAMMA_aux=(1./var_pi)-(1./Ve);
                gamma_aux=(nu_pi./var_pi)-(me./Ve);
                
                GAMMAupdated=beta*GAMMA_aux+(1-beta)*lambda;
                gammaupdated=beta*gamma_aux+(1-beta)*gamma;
                
                for k=1:window
                    if GAMMAupdated(k)>0 && sum(isfinite(cavity_i(k,:)))==M
                        meanLDPC(kk-ntaps+k)=gammaupdated(k)./GAMMAupdated(k);
                        varLDPC(kk-ntaps+k)=1./GAMMAupdated(k);
                    end
                end
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
    