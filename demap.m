function pvi=demap(pui,symbolmap,flaginterleaving,seed_inter)

% Function: demap
%
% pvi=demap(pui,symbolmap,flaginterleaving,seed_inter)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 27/02/2017
%
% Description: This function demaps and deinterleaves the probability of 
% symbols (pui) into the probability of bits (pvi)
% 
% Inputs: 
% pui is the probability of symbols
% symbolmap is a matrix where column i is the gray representation of symbol A(i)
% flaginterleaving is 1 if deinterleaving is included after demapping
% seed_inter is the random stream that determines the specific permutation
%
% Output: 
% pvi is the probability of bits

[L,N]=size(pui);
nbits=log2(L);
N=N*nbits;
pvi=zeros(2,N);

for k=1:nbits
    ind0=find(symbolmap(k,:)==0);
    ind1=find(symbolmap(k,:)==1);
    
    pvi(1,k:nbits:N)=sum(pui(ind0,:),1);
    pvi(2,k:nbits:N)=sum(pui(ind1,:),1);
end

% deinterleaving
if flaginterleaving
    pvi(1,:)=randdeintrlv(pvi(1,:),seed_inter);
    pvi(2,:)=randdeintrlv(pvi(2,:),seed_inter);
end
