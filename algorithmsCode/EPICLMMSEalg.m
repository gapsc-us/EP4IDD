function [x_decod,prob_b]=EPICLMMSEalg(A,complexFlag,sigma,y,H,pui_ldpc)

% Function: EPICLMMSEalg
%
% [x_decod,prob_b]=EPICLMMSEalg(A,complexFlag,sigma,y,H,pui_ldpc)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 28/06/2017
%
% Description: This function runs the EP IC-LMMSE [Senst11] for turbo 
% equalization. 
% It returns the estimated symbols (x_decod) and the extrinsic probability 
% for each symbol (prob_b)
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
% Outputs: 
% x_decod is the estimation of the transmitted symbols
% prob_b is the extrinsic probability for each symbol
%
% References: 
% [Senst11] = M. Senst and G. Ascheid, ?How the framework of expectation 
% propagation yields an iterative IC-LMMSE MIMO receiver?, in Proc. IEEE 
% Global Telecommunications Conference (GLOBECOM), Dec 2011, pp. 1-6. 

Nr=size(H,1);
Nx=size(H,2);
M=length(A);

AlphabetMat=(diag(A)*ones(M,Nx)).'; 

% Initialization
meq2dem=ones(Nx,M);
if isempty(pui_ldpc)
    pui_ldpc=1/M*ones(Nx,M);
end

%% First iteration
% The first iteration within the detector is equivalent to the conventional
% LMMSE

% Demapper Node Update
fx=meq2dem.*pui_ldpc; % approximated posterior
fx=fx./(diag(sum(fx,2))*ones(Nx,M)); % eq 33

[nufx,varfx]=meanvar(fx.',A);
if complexFlag
    % The complex system can be simulated by: 
    % - Considering the real and imaginary parts separately 
    % - Formulating the algorithms in complex formulation
    % The following is an adjustment to let the results of both
    % possibilities match: 
    sigma=sigma*sqrt(2);
    varfx=varfx.*2;
end

% This message is passed to the equalizer (eq 37)
mudem2eq=nufx;
vardem2eq=varfx;

% Equalizer node update
Cx=diag(vardem2eq);
Cy=H*Cx*H'+sigma^2*eye(Nr);
invCy=Cy\eye(Nr);

% Posterior
z=mudem2eq+Cx*H'*invCy*(y-H*mudem2eq); % eq 24 (or 26)
C=Cx-Cx*H'*invCy*H*Cx; % eq 25
varC=diag(C); % eq 27

% Extrinsic
vareq2dem=(varC.*vardem2eq)./(vardem2eq-varC); % eq 29
mueq2dem=(z.*vardem2eq-mudem2eq.*varC)./(vardem2eq-varC); % eq 28


%% Second iteration  
if complexFlag
    meq2dem=normpdfcomplex(AlphabetMat,mueq2dem*ones(1,M),sqrt(vareq2dem)*ones(1,M));
else
    meq2dem=normpdf(AlphabetMat,mueq2dem*ones(1,M),sqrt(vareq2dem)*ones(1,M));
end
meq2dem=meq2dem./(diag(sum(meq2dem,2))*ones(Nx,M));

% Demapper Node Update
fx=meq2dem.*pui_ldpc;
fx=fx./(diag(sum(fx,2))*ones(Nx,M)); % eq 33

[nufx,varfx]=meanvar(fx.',A);

% We now divide by meq2dem (eq 36)
mudem2eq=(nufx.*vareq2dem-mueq2dem.*varfx)./(vareq2dem-varfx);
vardem2eq=(varfx.*vareq2dem)./(vareq2dem-varfx);

% Looking for negative variances
indNeg=find(vardem2eq<0);
vardem2eq(indNeg)=varfx(indNeg);
mudem2eq(indNeg)=nufx(indNeg);

% Equalizer node update
Cx=diag(vardem2eq);
Cy=H*Cx*H'+sigma^2*eye(Nr);
invCy=Cy\eye(Nr);

% Posterior
z=mudem2eq+Cx*H'*invCy*(y-H*mudem2eq); % eq 24 (or 26)
C=Cx-Cx*H'*invCy*H*Cx; % eq 25
varC=diag(C); % eq 27

% Extrinsic
vareq2dem=(varC.*vardem2eq)./(vardem2eq-varC); % eq 29
mueq2dem=(z.*vardem2eq-mudem2eq.*varC)./(vardem2eq-varC); % eq 28

%% Final probability

if complexFlag
    prob_b=normpdfcomplex(AlphabetMat,mueq2dem*ones(1,M),sqrt(vareq2dem)*ones(1,M)).';
else
    prob_b=normpdf(AlphabetMat,mueq2dem*ones(1,M),sqrt(vareq2dem)*ones(1,M)).';
end

prob_b=prob_b./transpose(diag(sum(prob_b))*ones(Nx,M));
[~,ind]=max(prob_b);
x_decod=A(ind).';
