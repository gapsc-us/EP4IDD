clc,close all, clear all,
addpath('../')
addpath('../algorithmsCode')


run configFiles/MimoDEPconfig256QAM4x4Lim3Random

flagRun=1;

tic
mainEP(flagRun,dataEP,scenario)
toc

