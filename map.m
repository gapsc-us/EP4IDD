function pui=map(pvi,symbolmap,flaginterleaving,seed_inter)

% Function: map
%
% pui=map(pvi,symbolmap,flaginterleaving,seed_inter)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 16/01/2017
%
% Description: This function maps and interlaves the probability of bits 
% (pvi) into the probability of symbols (pui)
% 
% Inputs: 
% pvi is the probability of bits
% symbolmap is a matrix where column i is the gray representation of symbol A(i)
% flaginterleaving is 1 if interleaving is included after demapping
% seed_inter is the random stream that determines the specific permutation
%
% Output: 
% pui is the probability of symbols

[~,N]=size(pvi);
L=size(symbolmap,2);
nbits=log2(L);
Nu=N/nbits;
pui=ones(L,Nu);

% interleaving
if flaginterleaving
    pvi(1,:)=randintrlv(pvi(1,:),seed_inter);
    pvi(2,:)=randintrlv(pvi(2,:),seed_inter);
end

for kk=1:L
    indA=symbolmap(:,kk)+1;
    paux=ones(1,Nu);
    for k=1:nbits        
        paux=paux.*pvi(indA(k),k:nbits:N);
    end
    pui(kk,:)=paux;
end
