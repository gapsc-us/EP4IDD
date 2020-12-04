function [xe,LLR,varXe]=mmseTuchler(z, N1,N2,h,sigma,meanX,varX,energy)

% Function: mmseTuchler
%
% [xe,LLR,varXe]=mmseTuchler(z, N1,N2,h,sigma,meanX,varX,energy)
%
% Author: Juan José Murillo Fuentes
% Modified by Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 15/01/2018
%
% Description: This function runs the MMSE Wiener [Tuchler02a,Tuchler02b,
% Tuchler11] for turbo equalization. 
% It returns the estimated symbols (x_decod) and the extrinsic LLRs. This 
% function is used by DFEPalg.m, FMMSEalg.m and PFEPalg.m 
% 
% Inputs: 
% z is the received signal
% N2 is the number of symbols to be considered at the left-hand side of the
% window
% N1 is the number of symbols to be considered at the right-hand side of the
% window
% h is the channel
% sigma is the standard deviation of the noise
% meanX is the mean of the prior information 
% varX is the variance of the prior information 
% energy is the mean symbol energy of the constellation
%
% Output: 
% xe is the estimation of the transmitted symbols
% LLRs is the loglikehood ratio at the output of the equalizer
% varXe is the variance of the extrinsic distribution at the output of the 
% equalizer
%
% References: 
% equalization: principles and new results,? IEEE Trans. On Communications, 
% vol. 50, no. 5, pp. 754-767, May 2002.
% [Tuchler02b] = M. Tüchler, A. Singer, and R. Koetter, ?Minimum mean 
% squared error equalization using a priori information, ? IEEE Trans. On 
% Signal Processing, vol. 50, no 3, pp. 673-683, Mar 2002.
% [Tuchler11] = M. Tüchler, and A. Singer, ?Turbo equalization: An 
% overview,? IEEE Trans. On Information Theory, vol. 57, no. 2, pp. 920-952, 
% Feb 2011. 

h=flip(h);
M=length(h);
N=N1+N2+1;
Htc=obthTuchler(h,N);
s=Htc*[zeros(1,N2+M-1) 1 zeros(1,N1)]';

varXe=[];
Ns=length(z);
%Ampliamos todo metiendo zeros por izda y dcha
z=[zeros(N2+M,1); z; ones(N1,1)*sigma ]; %z=[randn(N2+M,1)*sigma; z; randn(N1,1)*sigma ];
meanX=[zeros(N2+M,1); meanX; zeros(N1,1)];
varX=[ones(N2+M,1)*1e-8; varX; ones(N1,1)*1e-8];
xe=[];
LLR=[];

for k1=1+M+N2:Ns+M+N2 %desplazo todo M+N2 a la derecha
    znv=z(k1-N2:k1+N1);
    mXv=meanX(k1-M-N2+1:k1+N1);
    Vn=diag(varX(k1-M-N2+1:k1+N1));
    cH=(sigma^2*eye(N)+Htc*Vn*Htc'+(energy-varX(k1))*(s*s'))\(s*energy); %
    xne=cH'*(znv-Htc*mXv+meanX(k1)*s);
    
    mun1=cH'*s;
    sigma2n1=cH'*s*energy*(1-s'*cH);
    LLR=[LLR; 2*mun1*xne/sigma2n1];
    xe=[xe; xne/mun1];
    varXe=[varXe; sigma2n1/(mun1^2)];
end
end

function H=obthTuchler(h,N)

% It obtains the matrix Hk of the filtered model yk=Hk*xk+wk, that will
% allow to compute the output xk of the channel 

h=h(:);
L=length(h);
H=zeros(N+L,N);

for i=1:N
    H(i,i:i+L-1)=h;
end

H=H(1:N,:);
end
