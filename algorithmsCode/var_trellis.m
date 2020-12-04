function [convmat_trans,convmat_salidas,transp_llegan_a_q,transq_llegan_a_p]=var_trellis(A,ntaps,heval)

% Function: var_trellis
%
% [convmat_trans,convmat_salidas,transp_llegan_a_q,transq_llegan_a_p]=var_trellis(A,ntaps,heval)
%
% Author: Luis Salamanca Miño
% Modified by Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 22/05/2016
%
% Description: This function obtains the trellis that will be used by
% BCJRalg.m to compute the BCJR algorithm

M=length(A); % constallation size

dec_states=0:M-1;
bin_states_aux=dec2bin(dec_states);
bin_states=zeros(M,log2(M));
for k1=1:M
    for k2=1:log2(M)
        bin_states(k1,k2)=str2num(bin_states_aux(k1,k2));
    end
end
all_states=zeros(M^(ntaps-1),ntaps-1);
all_bin_states=zeros(M^(ntaps-1),(ntaps-1)*log2(M));


% Matrix that keeps state transitions inside the trellis. Using this matrix
% we will be able to compute the other needed parameters. 
convmat_trans=zeros(M^(ntaps-1),M^(ntaps-1));

% Row-vector that will allow us to compute the previous matrix. It will
% be defined as the first row of the matrix, and then we apply a bit shift
% to the right per M rows. 
vect_rell=zeros(1,M^(ntaps-1));
vect_rell(1:M^(ntaps-2):M^(ntaps-1))=A;

% We fill the matrix
mat_aux_ones=ones(M^(ntaps-1),M);
for i=1:M:(M^(ntaps-1)-1)
    convmat_trans(i:i+M-1,:)=(diag(vect_rell)*mat_aux_ones).';
    vect_rell=[0 vect_rell(1:(M^(ntaps-1)-1))];
end

% We consider that p-states are referred to instant t and q-states to
% instant t+1 to which we evolve with a certain bit at the input of the
% channel

% For certain row-value, that indicates a q-state to which we move to, this
% matrix keeps the M-values of the p-states from which we can reach the
% q-state
transp_llegan_a_q=zeros(M^(ntaps-1),M);

% For certain row-value, that indicas a p-state to which we move to, this
% matrix keeps the M-values of the q-states from which we can reach the
% p-state
transq_llegan_a_p=zeros(M^(ntaps-1),M);

% convmat_salidas keeps, depending on the the current p-state, the value at
% the output withouth noise of an evolution to the q-state (that depends on
% the value of a bit at the input of the channel)
convmat_salidas=zeros(M^(ntaps-1),M^(ntaps-1));
sign_est=A(1)*ones(1,ntaps-1);

% Loop for computing the matrix convmat_salidas
for i=1:M^(ntaps-1) 
    % we fix the i-th column and keep the indeces of the rows with values
    % different to NaN
    [fila,~]=find(convmat_trans(:,i));
    transp_llegan_a_q(i,:)=fila;
    % we fix the i-th row and keep the indeces of the columns with values
    % different to NaN
    [~,columna]=find(convmat_trans(i,:));
    transq_llegan_a_p(i,:)=columna;
    
    % To sum up, what we do here is to obtain the chain of bits that
    % corresponds to the current state (k). By multiplying (element by
    % element) this chain and the heval from 2 to ntaps, we have the value
    % at the output of the channel. We have to include also the value
    % A(i)*heval(1), that depends on the element at the input of the
    % channel
    for kk=(ntaps-1):-1:1
        if (mod(i-1,M^(ntaps-1-kk))==0 && i~=1)
            ind=find(A==sign_est(kk)); % Index of alphabet element
            if ind==M
                ind=0;
            end
            sign_est(kk)=A(ind+1);
        end
        ind=find(A==sign_est(kk));
        all_bin_states(i,(kk-1)*log2(M)+1:(kk-1)*log2(M)+log2(M))=bin_states(ind,:);
    end
    all_states(i,:)=sign_est;
    for kk=1:M
        convmat_salidas(i,transq_llegan_a_p(i,kk))=A(kk)*heval(1,1)+sum(sign_est.*heval(2:ntaps,1).');
    end
end




