close all, clear all,
addpath('../')
addpath('../algorithmsCode')


run configFiles/MimoDEPconfig256QAM3x3Lim3Random

flagRun=1;

tic
mainEP(flagRun,dataEP,scenario)
toc

