function H=obth(h,N)

% Function: obth
%
% H=obth(h,N)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 22/05/2016
%
% Description: This function obtains the matrix H of the lineal model
% y=H*x+w, that will allow to compute the output x of the channel
% 
% Inputs: 
% h is a 1-by-L vector with the channel
% N is the length of the frame  
%
% Output: 
% H is a N+L-1-by-N matrix with the channel matrix

h=h(:);
L=length(h);
H=zeros(N+L,N);

for i=1:N
    H(i:i+L-1,i)=h;
end

H=H(1:N,1:N);