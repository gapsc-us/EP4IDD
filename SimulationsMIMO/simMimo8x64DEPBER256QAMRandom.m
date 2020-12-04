close all, clear all,
addpath('../')
addpath('../algorithmsCode')


run configFiles/MimoDEPconfig256QAM8x64Lim3Random

flagRun=1;

tic
mainEP(flagRun,dataEP,scenario)
toc

