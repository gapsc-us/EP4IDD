close all, clear all,
addpath('../')
addpath('../algorithmsCode')


run configFiles/MimoDEPconfig256QAM16x64Lim3Long25Random

flagRun=1;

tic
mainEP(flagRun,dataEP,scenario)
toc

