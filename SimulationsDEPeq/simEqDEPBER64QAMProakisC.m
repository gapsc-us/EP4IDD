close all, clear all,
addpath('../')
addpath('../algorithmsCode')


run configFiles/DEPconfig64QAMLim3ProakisC

flagRun=1;

tic
mainEP(flagRun,dataEP,scenario)
toc

