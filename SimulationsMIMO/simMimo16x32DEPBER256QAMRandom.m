close all, clear all,
addpath('../')
addpath('../algorithmsCode')


run configFiles/MimoDEPconfig256QAM16x32Lim3Random

flagRun=1;

tic
mainEP(flagRun,dataEP,scenario)
toc

