function [hard_detect_bcjr,prob_b_r]=BCJRalg(A,complexFlag,tam_bloque,h,x,n,sigma,convmat_trans,convmat_salidas,transp_llegan_a_q,transq_llegan_a_p,pui)

% Function: BCJRalg
%
% [hard_detect_bcjr,prob_b_r]=BCJRalg(A,complexFlag,tam_bloque,h,x,n,sigma,convmat_trans,convmat_salidas,transp_llegan_a_q,transq_llegan_a_p,pui)
%
% Author: Luis Salamanca Miño
% Modified by Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 22/05/2016
%
% Description: This function runs the BCJR [Bahl74] for turbo equalization. 
% It returns the estimated symbols (hard_detect_bcjr) and the probability 
% for each symbol (prob_b_r). If flagextrinsic=1, then prob_b is the 
% extrinsic distribution. If flagextrinsic=0, then prob_b is the posterior
% distribution. 
%
% Reference: 
% [Bahl74] = L. Bahl, J. Cocke, F. Jelinek, and J. Raviv, ?Optimal decoding 
% of linear codes for minimizing symbol error rate (corresp.),? IEEE Trans. 
% On Information Theory, vol. 20, no. 2, pp. 284-287, 1974. 

flagextrinsic=1; 

ntaps=length(h);
M=length(A);

% Obtención de la matrix para llevar a cabo la convolución con el canal
xmat=obthx(A(1),x.',length(h));

% Salida del canal
r=xmat*h.'+n;

% Para turbo estudio las probabilidades de simbolo
output_length=length(x);
input_length=output_length-ntaps+1;
if isempty(pui)
    pprior=[(1/M)*ones(M,input_length) [ones(1,ntaps-1);zeros(M-1,ntaps-1)]];
else
    pprior=[pui [ones(1,ntaps-1);zeros(M-1,ntaps-1)]];
end

%% FORWARD
% Initialization of matrices alpha, gamma_resul and prob_b_r

alphamat=zeros(M^(ntaps-1),tam_bloque);
gamma_resul=zeros(M^(ntaps-1)*M,tam_bloque);
prob_b_r=zeros(M,tam_bloque);

% We start by computing the alpha values as:
% alpha_t+1(q)=sum_p(alpha_t(p)*gamma_t(p,q))

alphamat(1,1,:)=1; % Consideramos que partimos del zero-state

for t=2:length(r) % From 2 because the first alpha is already computed 
    for p=1:M^(ntaps-1) % We compute all possible values of gamma_t(p,q)
        for i=1:M
            if complexFlag
                gamma_resul((p-1)*M+i,t-1)=pprior(i,t-1)*normpdfcomplex(r(t-1),convmat_salidas(p,transq_llegan_a_p(p,i)),sqrt(2)*sigma);
            else
                gamma_resul((p-1)*M+i,t-1)=pprior(i,t-1)*normpdf(r(t-1),convmat_salidas(p,transq_llegan_a_p(p,i)),sigma);
            end
        end
    end
    
    % Whis this matrix we know from which row we have to read gamma_resul,
    % depending on what was the previous state transition. In the rows we 
    % have "p" (current state) and, in the columns we have "q" (next
    % state). They match in the row's value of gamma_resul
    
    for q=1:M^(ntaps-1)
        for i=1:M
            alpha_i=alphamat(transp_llegan_a_q(q,i),t-1)...
                *gamma_resul((transp_llegan_a_q(q,i)-1)*M+find(A==convmat_trans(transp_llegan_a_q(q,i),q)),t-1);
            alphamat(q,t)=alphamat(q,t)+alpha_i;
        end
    end
    
    % Normalization of alpha in alphamat
    normalpha=sum(alphamat(:,t));
    alphamat(:,t)=alphamat(:,t)./normalpha;
end

%% BACKWARD
% Computing the value of gamma for the last time instant, which will be
% necessary in the backward procedure.

for p=1:M^(ntaps-1) % Computing all possible values for gamma_t(p,q)
    for i=1:M
        if complexFlag
            gamma_resul((p-1)*M+i,length(r))=pprior(i,length(r))*normpdfcomplex(r(length(r)),convmat_salidas(p,transq_llegan_a_p(p,i)),sqrt(2)*sigma);
        else
            gamma_resul((p-1)*M+i,length(r))=pprior(i,length(r))*normpdf(r(length(r)),convmat_salidas(p,transq_llegan_a_p(p,i)),sigma);
        end
    end
end

% Now we compute the beta values as:
% beta_t(p)=sum_q(beta_t+1(q)*gamma_t(p,q))

betaaux=zeros(M^(ntaps-1),1); % this vector keeps the values of beta at 
% instant t+1, with respect to the instant t (current instant in the loop)
betanew=zeros(M^(ntaps-1),1);
betaaux(1,1,:)=1; % We consider that the algorithm finishes in the zero-state
gamma_resul_aux=zeros(M^(ntaps-1)*M,1);

for t=length(r):-1:1
    for p=1:M^(ntaps-1)
        betanew(p,1)=0;
        for i=1:M
            beta_i=betaaux(transq_llegan_a_p(p,i),1)...
                *gamma_resul((p-1)*M+find(A==convmat_trans(p,transq_llegan_a_p(p,i))),t);
            betanew(p,1)=betanew(p,1)+beta_i;
        end
    end
    
    % We sort the state evolutions in matrix gamma_result, hence the first
    % M^(ntaps-1)/M positions correspond to evolutions because of a A(1) in
    % the input, the next ones to A(2) and so on until A(M)
    for entrada=1:M
        for p=1:M^(ntaps-1) % Now gamma_results keeps p(s_t=p,s_t+1=q|r)
            gamma_resul_aux((entrada-1)*M^(ntaps-1)+p)=alphamat(p,t)...
                *gamma_resul((p-1)*M+entrada,t)...
                *betaaux(transq_llegan_a_p(p,entrada),1);
        end
    end
    gamma_resul(:,t)=gamma_resul_aux;
    
    % Normalization of betanew, that will be used then to compute
    % everything
    betanorm=sum(betanew(:,1));
    betanew(:,1)=betanew(:,1)./betanorm;
    
    betaaux(:,1)=betanew(:,1);
    
    for entr=1:M
        prob_b_r(entr,t)=sum(gamma_resul((1+(entr-1)*M^(ntaps-1)):((M^(ntaps-1))+(entr-1)*M^(ntaps-1)),t));
    end
    
    prob_b_r_norm=sum(prob_b_r(:,t));
    prob_b_r(:,t)=prob_b_r(:,t)./prob_b_r_norm;
end

%% Probabilities at the output of the equalizer

[~,ind]=max(prob_b_r,[],1);
hard_detect_bcjr=A(ind).';

% Soft detection (extrinsics)
if flagextrinsic
    % Extrinsic probability
    prob_b_r_ext=prob_b_r(:,1:input_length)./pprior(:,1:input_length);
    K=sum(prob_b_r_ext);
    prob_b_r_ext=prob_b_r_ext./transpose(diag(K)*ones(input_length,M));
    prob_b_r=prob_b_r_ext;
end
end

function h=obthx(Ao,u,longeq)

output_length=length(u);

h=zeros(output_length,longeq);

for i=1:longeq;
    h(:,i)=[Ao*ones(1,i-1) u(1,1:(output_length-i+1))].';
end;
end