function [x_decod,prob_b,x_decod_FB,prob_b_FB]=SEPalg(A,complexFlag,sigma,r,h,pui_ldpc)

% Function: SEPalg
%
% [x_decod,prob_b,x_decod_FB,prob_b_FB]=SEPalg(A,complexFlag,sigma,r,h,pui_ldpc)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 8/02/2018
%
% Description: This function runs the SEP [Santos17] for turbo equalization. 
% It returns the estimated symbols (x_decod) and the probability for each 
% symbol (prob_b). If flagextrinsic=1, then prob_b is the extrinsic 
% distribution. If flagextrinsic=0, then prob_b is the posterior
% distribution. 
% It also computes the estimation (x_decod_FB) and probabilities (prob_b_FB)
% in the backward and forward steps. 
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
%
% Output: 
% x_decod is the estimation of the transmitted symbols
% prob_b is the (posterior or extrinsic) probability for each symbol after
% the smoothing step
% x_decod_FB is the estimation of the transmitted symbols after the
% forward/backward steps
% prob_b_FB is the (posterior or extrinsic) probability for each symbol after
% the forward/backward steps
%
% References: 
% [Santos17] = I. Santos, J. J. Murillo-Fuentes, E. Arias-de-Reyna, and P. 
% M. Olmos, ?Probabilistic equalization with a smoothing expectation 
% propagation approach,? IEEE Transactions on Wireless Communications, vol. 
% 16, no. 5, pp. 2950-2962, May 2017. 

flagextrinsic=1; % to 1 if the probability sent to the decoder is the extrinsic (cavity). to 0 if it is the posterior. 

M=length(A);
ntaps=length(h);
output_length=length(r);
input_length=output_length-ntaps+1;
AlphabetMat=(diag(A)*ones(M,input_length)).'; 

varmin=1e-3;
A0=0;

%% EP parameters in [Santos17]
if isempty(pui_ldpc)   
    L_FBEP=10;
else
    L_FBEP=2;
end
L_EP=0;
eps_mep=0.5;
eps_soft=1e-8;
beta=0.1;

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
    varLDPC=varLDPC*2;
    sigma=sigma*sqrt(2);
end

pui=[meanLDPC varLDPC];

Hall=obth(h,2*ntaps);
Hmat=Hall(ntaps:end-1,1:end-1);

h_flip=fliplr(h);

x_decod=zeros(input_length,1);

% Backward
x_decod_b=zeros(input_length,1);

gamma=(A0/varmin)*ones(2*ntaps-1,1); % Modif 21 feb 17
lambda=(1/varmin)*ones(2*ntaps-1,1);
C_q=varmin*eye(2*ntaps-1);
nu_q=A0*ones(2*ntaps-1,1); % Modif 21 feb 17

% Mean and variance of q in each iteration
nu_qback=zeros(2*ntaps-1,output_length);
Cinv_qback=zeros(2*ntaps-1,2*ntaps-1,output_length);
% Gamma and lambda of each iteration in backward steps
gamma_b=zeros(2*ntaps-1,output_length);
lambda_b=zeros(2*ntaps-1,output_length);

%% Backwards

nub_b=zeros(input_length,1);
varb_b=zeros(input_length,1);

% Hacemos el backwards como un forward dándole la vuelta a los taps, u e y.
hb=fliplr(h);
hb_flip=fliplr(hb);
rb=flipud(r);
puib=flipud(pui);
if ~isempty(pui)    
    puibaux=[A0*ones(2*ntaps-1,1) varmin*ones(2*ntaps-1,1)];
    puiaux=[A0*ones(2*ntaps-1,1) varmin*ones(2*ntaps-1,1)];
else
    puibaux=[];
    puiaux=[];
end

for kk=1:output_length
    if kk<2*ntaps-1 
        if ~isempty(pui)
            puibaux=[puibaux(2:end,:);puib(kk,:)];
        end

        Cq_inv=inv(C_q(2:end,2:end));
        H=hb_flip;
        C_inv=sigma^(-2)*[zeros(ntaps-1,2*ntaps-1);zeros(ntaps,ntaps-1) H'*H]+[Cq_inv-diag(lambda(2:end)) zeros(2*ntaps-2,1); zeros(1,2*ntaps-1)];
        zCinv=sigma^(-2)*[zeros(ntaps-1,1); H'*rb(kk)]+[Cq_inv*nu_q(2:end)-gamma(2:end);0];
        
        [xaux,gamma,lambda,nu_q,C_q]=EP_LD(A,complexFlag,2*ntaps-1,L_FBEP,zCinv,C_inv,gamma(2:end),lambda(2:end),eps_mep,beta,2*ntaps-1-kk+1:2*ntaps-1,puibaux);
        x_decod_b(input_length-kk+1:input_length)=flipud(xaux(2*ntaps-1-kk+1:2*ntaps-1));
        nub_b(input_length-kk+1:input_length)=flipud(nu_q(2*ntaps-1-kk+1:2*ntaps-1));
        varb_b(input_length-kk+1:input_length)=flipud(diag(C_q(2*ntaps-1-kk+1:2*ntaps-1,2*ntaps-1-kk+1:2*ntaps-1)));
              
        % BCJR
        nu_qback(:,output_length-kk+1)=nu_q(2*ntaps-1:-1:1);
        Cinv_qback(:,:,output_length-kk+1)=inv(C_q(2*ntaps-1:-1:1,2*ntaps-1:-1:1));
        
        gamma_b(:,output_length-kk+1)=gamma(2*ntaps-1:-1:1);
        lambda_b(:,output_length-kk+1)=lambda(2*ntaps-1:-1:1);
            
    elseif kk<input_length+1
        if ~isempty(pui)
            puibaux=[puibaux(2:end,:);puib(kk,:)];
        end

        Cq_inv=inv(C_q(2:end,2:end));
        H=hb_flip;
        C_inv=sigma^(-2)*[zeros(ntaps-1,2*ntaps-1);zeros(ntaps,ntaps-1) H'*H]+[Cq_inv-diag(lambda(2:end)) zeros(2*ntaps-2,1); zeros(1,2*ntaps-1)];
        zCinv=sigma^(-2)*[zeros(ntaps-1,1); H'*rb(kk)]+[Cq_inv*nu_q(2:end)-gamma(2:end);0];

        [x_decod_b(output_length-kk+1+ntaps-1:-1:input_length-kk+1),gamma,lambda,nu_q,C_q]=EP_LD(A,complexFlag,2*ntaps-1,L_FBEP,zCinv,C_inv,gamma(2:end),lambda(2:end),eps_mep,beta,1:2*ntaps-1,puibaux);
        nub_b(input_length-kk+1:output_length-kk+1+ntaps-1)=flipud(nu_q);
        varb_b(input_length-kk+1:output_length-kk+1+ntaps-1)=flipud(diag(C_q));

        % BCJR
        nu_qback(:,output_length-kk+1)=nu_q(2*ntaps-1:-1:1);
        Cinv_qback(:,:,output_length-kk+1)=inv(C_q(2*ntaps-1:-1:1,2*ntaps-1:-1:1));
        
        gamma_b(:,output_length-kk+1)=gamma(2*ntaps-1:-1:1);
        lambda_b(:,output_length-kk+1)=lambda(2*ntaps-1:-1:1);
    else
        if ~isempty(pui)
            % Modif 21 feb 17
            puibaux=[puibaux(2:end,:); [A0 varmin]];
        end
        
        Cq_inv=inv(C_q(2:end,2:end));
        H=hb_flip;
        C_inv=sigma^(-2)*[zeros(ntaps-1,2*ntaps-1);zeros(ntaps,ntaps-1) H'*H]+[Cq_inv-diag(lambda(2:end)) zeros(2*ntaps-2,1); zeros(1,2*ntaps-1)];
        zCinv=sigma^(-2)*[zeros(ntaps-1,1); H'*rb(kk)]+[Cq_inv*nu_q(2:end)-gamma(2:end);0];

        % Modif 21 feb 2017
        [xaux,gamma,lambda,nu_q,C_q]=EP_LD(A,complexFlag,2*ntaps-1,L_FBEP,zCinv,C_inv,[gamma(2:end);A0/varmin],[lambda(2:end);1/varmin],eps_mep,beta,1:2*ntaps-1-kk+input_length,puibaux);
        x_decod_b(output_length-kk+1+ntaps-1:-1:1)=xaux(1:2*ntaps-1-kk+input_length);
        nub_b(output_length-kk+1+ntaps-1:-1:1)=nu_q(1:2*ntaps-1-kk+input_length);
        varb_b(output_length-kk+1+ntaps-1:-1:1)=diag(C_q(1:2*ntaps-1-kk+input_length,1:2*ntaps-1-kk+input_length));
         
        % BCJR
        nu_qback(:,output_length-kk+1)=nu_q(2*ntaps-1:-1:1);
        Cinv_qback(:,:,output_length-kk+1)=inv(C_q(2*ntaps-1:-1:1,2*ntaps-1:-1:1));
        
        gamma_b(:,output_length-kk+1)=gamma(2*ntaps-1:-1:1);
        lambda_b(:,output_length-kk+1)=lambda(2*ntaps-1:-1:1);
    end
end


%% Forward and smoothing

x_decod_f=zeros(input_length,1);
nub_f=zeros(input_length,1);
varb_f=zeros(input_length,1);

nub=zeros(input_length,1);
varb=zeros(input_length,1);

gamma=(A0/varmin)*ones(2*ntaps-1,1); 
lambda=(1/varmin)*ones(2*ntaps-1,1);
C_q=varmin*eye(2*ntaps-1);
nu_q=A0*ones(2*ntaps-1,1); 

for kk=1:output_length
    if kk<2*ntaps-1
        if ~isempty(pui)
            puiaux=[puiaux(2:end,:);pui(kk,:)];
        end
        
        Cq_inv=inv(C_q(2:end,2:end));
        H=h_flip;
        C_inv=sigma^(-2)*[zeros(ntaps-1,2*ntaps-1);zeros(ntaps,ntaps-1) H'*H]+[Cq_inv-diag(lambda(2:end)) zeros(2*ntaps-2,1); zeros(1,2*ntaps-1)];
        zCinv=sigma^(-2)*[zeros(ntaps-1,1); H'*r(kk)]+[Cq_inv*nu_q(2:end)-gamma(2:end);0];
        
        [xaux,gamma,lambda,nu_q,C_q]=EP_LD(A,complexFlag,2*ntaps-1,L_FBEP,zCinv,C_inv,gamma(2:end),lambda(2:end),eps_mep,beta,2*ntaps-1-kk+1:2*ntaps-1,puiaux);
        x_decod_f(1:kk)=xaux(2*ntaps-1-kk+1:2*ntaps-1);
        nub_f(1:kk)=nu_q(2*ntaps-1-kk+1:2*ntaps-1);
        varb_f(1:kk)=diag(C_q(2*ntaps-1-kk+1:2*ntaps-1,2*ntaps-1-kk+1:2*ntaps-1));
                
        if kk>=ntaps
            % Denominator
            H=Hmat;
            Cinv_p=sigma^(-2)*H'*H;
            zCinv_p=sigma^(-2)*H'*r(kk-ntaps+1:kk);
            gamma_den=[gamma_b(1:ntaps,kk-ntaps+1);gamma(ntaps+1:end)];
            lambda_den=[lambda_b(1:ntaps,kk-ntaps+1);lambda(ntaps+1:end)];
            Covinv_den=Cinv_p+diag(lambda_den);
            Covinvnu_den=zCinv_p+gamma_den;
            
            % Numerator
            C_qinvf=inv(C_q);
            Covinv_num=C_qinvf+Cinv_qback(:,:,kk-ntaps+1);
            Covinvnu_num=C_qinvf*nu_q+Cinv_qback(:,:,kk-ntaps+1)*nu_qback(:,kk-ntaps+1);
            
            % FWEP
            Cinv_qall=Covinv_num-Covinv_den;
            Cinvnu_qall=Covinvnu_num-Covinvnu_den;
            
            gamma_all=[gamma(1:ntaps);gamma_b(ntaps+1:end,kk-ntaps+1)];
            lambda_all=[lambda(1:ntaps);lambda_b(ntaps+1:end,kk-ntaps+1)];
            % Dividimos entre el producto de exponenciales
            Cinv_fb_exp=Cinv_qall-diag(lambda_all);
            Cnu_fb_exp=Cinvnu_qall-gamma_all;
            
            [xtotal,~,~,nutotal,Ctotal,me,Ve]=EP_LD(A,complexFlag,2*ntaps-1,L_EP,Cnu_fb_exp,Cinv_fb_exp,gamma_all,lambda_all,eps_soft,beta,2*ntaps-kk:2*ntaps-1,[]);
            
            x_decod(kk-ntaps+1)=xtotal(ntaps);
            if flagextrinsic
                nub(kk-ntaps+1)=me(ntaps);
                varb(kk-ntaps+1)=Ve(ntaps);
            else
                nub(kk-ntaps+1)=nutotal(ntaps);
                varb(kk-ntaps+1)=Ctotal(ntaps,ntaps);
            end
        end
    elseif kk<=input_length
        if ~isempty(pui)
            puiaux=[puiaux(2:end,:);pui(kk,:)];
        end
        
        Cq_inv=inv(C_q(2:end,2:end));
        H=h_flip;
        C_inv=sigma^(-2)*[zeros(ntaps-1,2*ntaps-1);zeros(ntaps,ntaps-1) H'*H]+[Cq_inv-diag(lambda(2:end)) zeros(2*ntaps-2,1); zeros(1,2*ntaps-1)];
        zCinv=sigma^(-2)*[zeros(ntaps-1,1); H'*r(kk)]+[Cq_inv*nu_q(2:end)-gamma(2:end);0];
        
        [x_decod_f(kk-2*ntaps+2:kk),gamma,lambda,nu_q,C_q]=EP_LD(A,complexFlag,2*ntaps-1,L_FBEP,zCinv,C_inv,gamma(2:end),lambda(2:end),eps_mep,beta,1:2*ntaps-1,puiaux);
        nub_f(kk-2*ntaps+2:kk)=nu_q;
        varb_f(kk-2*ntaps+2:kk)=diag(C_q);
        
        % Denominator
        H=Hmat;
        Cinv_p=sigma^(-2)*H'*H;
        zCinv_p=sigma^(-2)*H'*r(kk-ntaps+1:kk);
        gamma_den=[gamma_b(1:ntaps,kk-ntaps+1);gamma(ntaps+1:end)];
        lambda_den=[lambda_b(1:ntaps,kk-ntaps+1);lambda(ntaps+1:end)];
        Covinv_den=Cinv_p+diag(lambda_den);
        Covinvnu_den=zCinv_p+gamma_den;
        
        % Numerator
        C_qinvf=inv(C_q);
        Covinv_num=C_qinvf+Cinv_qback(:,:,kk-ntaps+1);
        Covinvnu_num=C_qinvf*nu_q+Cinv_qback(:,:,kk-ntaps+1)*nu_qback(:,kk-ntaps+1);
        
        % FWEP
        Cinv_qall=Covinv_num-Covinv_den;
        Cinvnu_qall=Covinvnu_num-Covinvnu_den;
        gamma_all=[gamma(1:ntaps);gamma_b(ntaps+1:end,kk-ntaps+1)];
        lambda_all=[lambda(1:ntaps);lambda_b(ntaps+1:end,kk-ntaps+1)];
        % Dividimos entre el producto de exponenciales
        Cinv_fb_exp=Cinv_qall-diag(lambda_all);
        Cnu_fb_exp=Cinvnu_qall-gamma_all;
        
        [xtotal,~,~,nutotal,Ctotal,me,Ve]=EP_LD(A,complexFlag,2*ntaps-1,L_EP,Cnu_fb_exp,Cinv_fb_exp,gamma_all,lambda_all,eps_soft,beta,1:2*ntaps-1,[]);
        
        x_decod(kk-ntaps+1)=xtotal(ntaps);
        if flagextrinsic
            nub(kk-ntaps+1)=me(ntaps);
            varb(kk-ntaps+1)=Ve(ntaps);
        else
            nub(kk-ntaps+1)=nutotal(ntaps);
            varb(kk-ntaps+1)=Ctotal(ntaps,ntaps);
        end
    else
        if ~isempty(pui)
            puiaux=[puiaux(2:end,:);[A0 varmin]];
        end
        
        Cq_inv=inv(C_q(2:end,2:end));
        H=h_flip;
        C_inv=sigma^(-2)*[zeros(ntaps-1,2*ntaps-1);zeros(ntaps,ntaps-1) H'*H]+[Cq_inv-diag(lambda(2:end)) zeros(2*ntaps-2,1); zeros(1,2*ntaps-1)];
        zCinv=sigma^(-2)*[zeros(ntaps-1,1); H'*r(kk)]+[Cq_inv*nu_q(2:end)-gamma(2:end);0];
        
        [xaux,gamma,lambda,nu_q,C_q]=EP_LD(A,complexFlag,2*ntaps-1,L_FBEP,zCinv,C_inv,gamma(2:end),lambda(2:end),eps_mep,beta,1:2*ntaps-1-kk+input_length,puiaux);
        x_decod_f(kk-2*ntaps+2:input_length)=xaux(1:2*ntaps-1-kk+input_length);
        nub_f(kk-2*ntaps+2:input_length)=nu_q(1:2*ntaps-1-kk+input_length);
        varb_f(kk-2*ntaps+2:input_length)=diag(C_q(1:2*ntaps-1-kk+input_length,1:2*ntaps-1-kk+input_length));

        % Denominator
        H=Hmat;
        Cinv_p=sigma^(-2)*H'*H;
        zCinv_p=sigma^(-2)*H'*r(kk-ntaps+1:kk);
        gamma_den=[gamma_b(1:ntaps,kk-ntaps+1);gamma(ntaps+1:end)];
        lambda_den=[lambda_b(1:ntaps,kk-ntaps+1);lambda(ntaps+1:end)];
        Covinv_den=Cinv_p+diag(lambda_den);
        Covinvnu_den=zCinv_p+gamma_den;
        
        % Numerator
        C_qinvf=inv(C_q);
        Covinv_num=C_qinvf+Cinv_qback(:,:,kk-ntaps+1);
        Covinvnu_num=C_qinvf*nu_q+Cinv_qback(:,:,kk-ntaps+1)*nu_qback(:,kk-ntaps+1);
        
        % FWEP
        Cinv_qall=Covinv_num-Covinv_den;
        Cinvnu_qall=Covinvnu_num-Covinvnu_den;
        gamma_all=[gamma(1:ntaps);gamma_b(ntaps+1:end,kk-ntaps+1)];
        lambda_all=[lambda(1:ntaps);lambda_b(ntaps+1:end,kk-ntaps+1)];
        % Dividimos entre el producto de exponenciales
        Cinv_fb_exp=Cinv_qall-diag(lambda_all);
        Cnu_fb_exp=Cinvnu_qall-gamma_all;

        [xtotal,~,~,nutotal,Ctotal,me,Ve]=EP_LD(A,complexFlag,2*ntaps-1,L_EP,Cnu_fb_exp,Cinv_fb_exp,gamma_all,lambda_all,eps_soft,beta,1:2*ntaps-1-kk+input_length,[]);
        
        x_decod(kk-ntaps+1)=xtotal(ntaps);
        if flagextrinsic
            nub(kk-ntaps+1)=me(ntaps);
            varb(kk-ntaps+1)=Ve(ntaps);
        else
            nub(kk-ntaps+1)=nutotal(ntaps);
            varb(kk-ntaps+1)=Ctotal(ntaps,ntaps);
        end        
    end
end

%% Probabilities at the output of the equalizer

% BCJR
if complexFlag
    prob_b=normpdfcomplex(AlphabetMat,nub*ones(1,M),sqrt(varb)*ones(1,M)).';
else
    prob_b=normpdf(AlphabetMat,nub*ones(1,M),sqrt(varb)*ones(1,M)).';
end
prob_b=prob_b./transpose(diag(sum(prob_b))*ones(input_length,M));

prob_b_FB=zeros(M,input_length,2);
x_decod_FB=zeros(input_length,2);
% FWEP
if complexFlag
    prob_b_FB(:,:,1)=normpdfcomplex(AlphabetMat,nub_f*ones(1,M),sqrt(varb_f)*ones(1,M)).';
else
    prob_b_FB(:,:,1)=normpdf(AlphabetMat,nub_f*ones(1,M),sqrt(varb_f)*ones(1,M)).';
end
K=sum(prob_b_FB(:,:,1));
prob_b_FB(:,:,1)=prob_b_FB(:,:,1)./transpose(diag(K)*ones(input_length,M));
x_decod_FB(:,1)=x_decod_f;
% BWEP
if complexFlag
    prob_b_FB(:,:,2)=normpdfcomplex(AlphabetMat,nub_b*ones(1,M),sqrt(varb_b)*ones(1,M)).';
else
    prob_b_FB(:,:,2)=normpdf(AlphabetMat,nub_b*ones(1,M),sqrt(varb_b)*ones(1,M)).';
end
K=sum(prob_b_FB(:,:,2));
prob_b_FB(:,:,2)=prob_b_FB(:,:,2)./transpose(diag(K)*ones(input_length,M));
x_decod_FB(:,2)=x_decod_b;
    
