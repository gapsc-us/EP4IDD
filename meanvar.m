function [pui_mean,pui_var]=meanvar(pui,A)

% Function: meanvar
%
% [pui_mean,pui_var]=meanvar(pui,A)
%
% Author: Irene Santos Velázquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 27/12/2016
%
% Description: This function computes the mean and variance of a
% probability of symbols (pui) 
% 
% Inputs: 
% pui is the probability of symbols
% A is the set of symbols 
%
% Outputs: 
% pui_mean is the mean of pui
% pui_var is the variance of pui

[M,input_length]=size(pui);
AlphabetMat=(diag(A)*ones(M,input_length)).'; 

pui_mean=sum((diag(A)*pui).',2);
pui_var=sum(conj(AlphabetMat-diag(pui_mean)*ones(input_length,M)).*(AlphabetMat-diag(pui_mean)*ones(input_length,M)).*pui',2);

% To avoid numerical instabilities we set a minimum variance of 1e-8
pui_var=max([pui_var,1e-8*ones(input_length,1)],[],2);