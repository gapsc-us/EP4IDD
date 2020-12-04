function [x_decod,prob_b]=ML(A,complexFlag,sigma,y,H,pui_ldpc)

% Function: ML
%
% [x_decod,prob_b]=ML(A,complexFlag,sigma,y,H,pui_ldpc)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 7/03/2017
%
% Description: This function runs the MAP for turbo equalization. 
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

flagextrinsic=1;

input_length=size(H,2);
M=length(A);

if isempty(pui_ldpc)
    pui_ldpc=1/M*ones(M,input_length);
end

% Decoder
prob_join=zeros(M^input_length,1);

% First word to be checked
indexword=ones(input_length,1);
x_ml=A(indexword).';
puiword=pui_ldpc(1,:);

% Initialitation of ML
min_prob=mvnpdf(y.',(H*x_ml).',sigma.^2*eye(input_length))*prod(puiword);
x_decod=x_ml;

for k=1:M^input_length
    for kk=input_length:-1:1
        if (mod(k-1,M^(input_length-kk))==0 && k~=1)
            ind=find(A==x_ml(kk)); % Index of alphabet element
            if ind==M
                ind=0;
            end
            indexword(kk)=ind+1;
            x_ml=A(indexword).';
            puiword(kk)=pui_ldpc(ind+1,kk);
        end
    end
    
    if complexFlag
        prob_join(k)=real(mvnpdfcomplex(y,(H*x_ml),2*sigma.^2*eye(input_length))*prod(puiword));
    else
        prob_join(k)=real(mvnpdf(y.',(H*x_ml).',sigma.^2*eye(input_length))*prod(puiword));
    end
    
    % ML solution
    prob_act=prob_join(k);
    if prob_act>min_prob
        x_decod=x_ml;
        min_prob=prob_act;
    end
end

prob_join=prob_join/sum(prob_join);

prob_b=pjoin2bit(prob_join,input_length,M);

if flagextrinsic
    % Extrinsic probability
    prob_b=prob_b./pui_ldpc;
    prob_b=prob_b./transpose(diag(sum(prob_b))*ones(input_length,M));
end
end

function pbit=pjoin2bit(pjoin,input_length,M)

pbit=zeros(M,input_length);
njoin=length(pjoin);

for col=1:input_length
    if col==1
        step=njoin/M;
    else
        step=step/M;
    end
    for row=1:M
        ind_ini=(row-1)*step+1;
        ind_fin=ind_ini+step-1;
        aux=0;
        while ind_ini<=njoin
            aux=aux+sum(pjoin(ind_ini:ind_fin,:),1);
            ind_ini=ind_ini+M*step;
            ind_fin=ind_fin+M*step;
        end   
        pbit(row,col)=aux;
    end    
end
end