close all, clear all,
addpath('../')
addpath('../algorithmsCode')

run configFiles/EPconfig256QAMLim3ProakisC
%run configFiles/EPconfig256QAMLim5ProakisC

flagRun=1;

tic
mainEP(flagRun,dataEP,scenario)
toc

