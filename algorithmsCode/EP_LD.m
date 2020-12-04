function [x_decod,gamma,GAMMA,nu_qi,C_q,ti,hi_2]=EP_LD(A,complexFlag,N,L_EP,zCinv,Cinv,gamma0,GAMMA0,eps0,beta,indReest,pui)

% Function: EP_LD
%
% [x_decod,gamma,GAMMA,nu_qi,C_q,ti,hi_2]=EP_LD(A,complexFlag,N,L_EP,zCinv,Cinv,gamma0,GAMMA0,eps0,beta,indReest,pui)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 8/02/2017
%
% Description: This function is used by SEPalg and computes the EP 
% algorithm applied in [Santos17]
%
% Reference: 
% [Santos17] = I. Santos, J. J. Murillo-Fuentes, E. Arias-de-Reyna, and P. 
% M. Olmos, ?Probabilistic equalization with a smoothing expectation 
% propagation approach,? IEEE Transactions on Wireless Communications, vol. 
% 16, no. 5, pp. 2950-2962, May 2017. 

M=length(A);
energy=sum(abs(A).^2)/M;
if complexFlag
    % The complex system can be simulated by: 
    % - Considering the real and imaginary parts separately 
    % - Formulating the algorithms in complex formulation
    % The following is an adjustment to let the results of both
    % possibilities match: 
    energy=energy*2;
end

AlphabetMat=(diag(A)*ones(M,N)).'; 

% Initialization
gamma=zeros(N,1);
GAMMA=zeros(N,1);
if isempty(gamma0)
    if isempty(pui)
        gamma=zeros(N,1);
        GAMMA=ones(N,1)/energy;
    else
        GAMMA=1./pui(:,2);
        gamma=pui(:,1).*GAMMA;
    end
elseif length(zCinv)==length(gamma0)+1
    GAMMA(1:N-1)=GAMMA0;
    gamma(1:N-1)=gamma0;
    if isempty(pui)  
        GAMMA(N)=1/energy;
    else
        GAMMA(N)=1./pui(N,2);
        gamma(N)=pui(N,1).*GAMMA(N);
    end
    
else
    gamma(1:N)=gamma0;
    GAMMA(1:N)=GAMMA0;
end


% If eps0<0 execute optimal EP. In eps0>0 execute fix EP
if eps0>0
    epsilon=eps0;
end

% meand and variance vectors of all marginal qi(x) at iteration
% l=0
C_q=inv(Cinv+diag(GAMMA));
var_qi=diag(C_q);
nu_qi=C_q*(zCinv+gamma);

% Mean (ti) and variance (hi_2) of each cavity marginal
hi_2=var_qi./(1-var_qi.*GAMMA);
ti=hi_2.*(nu_qi./var_qi-gamma);

for l=1:L_EP
    
    if eps0<0
        epsilon=max([2.^(-max([l-5, 1])),1e-8]);
    end 
    
    % Compute cavity marginal for each xi
    if complexFlag
        cavity_i=normpdfcomplex(AlphabetMat,ti*ones(1,M),sqrt(hi_2)*ones(1,M));
    else
        cavity_i=normpdf(AlphabetMat,ti*ones(1,M),sqrt(hi_2)*ones(1,M));
    end
    
    % Normalize distribution (33)
    K=sum(cavity_i,2);
    cavity_i=cavity_i./(diag(K)*ones(N,M));
    
    % Mean (nu_pi) and variance (var_pi) of distribution (33)
    nu_pi=sum((diag(A)*transpose(cavity_i)).',2);
    var_pi=sum(conj(AlphabetMat-diag(nu_pi)*ones(N,M)).*(AlphabetMat-diag(nu_pi)*ones(N,M)).*cavity_i,2);
    var_pi=max([var_pi,epsilon*ones(N,1)],[],2);
    
    % Compute new values for gamma and GAMMA
    GAMMA_aux=(1./var_pi)-(1./hi_2);
    gamma_aux=(nu_pi./var_pi)-(ti./hi_2);
    
    GAMMAupdated=beta*GAMMA_aux+(1-beta)*GAMMA;
    gammaupdated=beta*gamma_aux+(1-beta)*gamma;
    
    for k=1:length(indReest)
        if GAMMAupdated(indReest(k))>0 && sum(isfinite(cavity_i(indReest(k),:)))==M 
            GAMMA(indReest(k))=GAMMAupdated(indReest(k));
            gamma(indReest(k))=gammaupdated(indReest(k));
        end
    end
    
    % meand and variance vectors of all marginal qi(x) at next
    % iteration l
    C_q=inv(Cinv+diag(GAMMA));
    var_qi=diag(C_q);
    nu_qi=C_q*(zCinv+gamma);
    
    % Mean (ti) and variance (hi_2) of each cavity marginal
    hi_2=var_qi./(1-var_qi.*GAMMA);
    ti=hi_2.*(nu_qi./var_qi-gamma);
end

% Decoder
if complexFlag
    prob_b=normpdfcomplex(AlphabetMat,nu_qi*ones(1,M),sqrt(var_qi)*ones(1,M)).';
else
    prob_b=normpdf(AlphabetMat,nu_qi*ones(1,M),sqrt(var_qi)*ones(1,M)).';
end

prob_b=prob_b./transpose(diag(sum(prob_b))*ones(N,M));
[~,ind]=max(prob_b);
x_decod=A(ind).';