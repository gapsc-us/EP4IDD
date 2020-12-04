close all, clear all,
addpath('../')
addpath('../algorithmsCode')

%run configFiles/EPconfig256QAMLim5Random
run configFiles/EPconfig256QAMLim3Random

flagRun=1;

tic
mainEP(flagRun,dataEP,scenario)
toc

