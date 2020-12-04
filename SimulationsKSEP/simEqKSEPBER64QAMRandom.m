close all, clear all,
addpath('../')
addpath('../algorithmsCode')

%run configFiles/EPconfig64QAMLim5Random
run configFiles/EPconfig64QAMLim3Random

flagRun=1;

tic
mainEP(flagRun,dataEP,scenario)
toc

