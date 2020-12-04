close all, clear all,
addpath('../')
addpath('../algorithmsCode')


run configFiles/MimoDEPconfig64QAM16x64Lim3Random

flagRun=1;

tic
mainEP(flagRun,dataEP,scenario)
toc

