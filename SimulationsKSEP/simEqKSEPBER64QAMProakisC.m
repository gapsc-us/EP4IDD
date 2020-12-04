close all, clear all,
addpath('../')
addpath('../algorithmsCode')

%run configFiles/EPconfig64QAMLim5ProakisC.m
run configFiles/EPconfig64QAMLim3ProakisC.m
flagRun=1;

tic
mainEP(flagRun,dataEP,scenario)
toc
