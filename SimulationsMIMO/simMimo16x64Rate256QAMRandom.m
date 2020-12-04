close all,clear all
addpath('../')
addpath('../algorithmsCode')
%addpath('../ParityCheckMatrix')
%Used to generate and compare curves with different LLRlim and EPiterations
tic
flagRun=1;

%1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 5/6, 8/9, or 9/10.
dataEP.rate=1/4;
dataEP.SNRdBIni=4-3; %11.5-3;
dataEP.SNRdBEnd=7-3; %15-3
run configFiles/MimoDEP4Rateconfig256QAM16x64Lim3Random
disp('R=1/4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
mainEP(flagRun,dataEP,scenario),toc


dataEP.rate=1/3;
dataEP.SNRdBIni=6-3; %11.5-3;
dataEP.SNRdBEnd=10-3; %15-3
run configFiles/MimoDEP4Rateconfig256QAM16x64Lim3Random
disp('N=1/3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
mainEP(flagRun,dataEP,scenario),toc


dataEP.rate=1/2;
dataEP.SNRdBIni=11-3; %11.5-3;
dataEP.SNRdBEnd=15-3; %15-3
run configFiles/MimoDEP4Rateconfig256QAM16x64Lim3Random
disp('N=1/2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
mainEP(flagRun,dataEP,scenario),toc

dataEP.rate=3/4;
dataEP.SNRdBIni=17-3; %11.5-3;
dataEP.SNRdBEnd=20-3; %15-3
run configFiles/MimoDEP4Rateconfig256QAM16x64Lim3Random
disp('N=3/4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
mainEP(flagRun,dataEP,scenario),toc

% dataEP.rate=4/5;
% dataEP.SNRdBIni=17-3; %11.5-3;
% dataEP.SNRdBEnd=20-3; %15-3
% run configFiles/MimoDEP4Rateconfig256QAM16x64Lim3Random
% disp('N=4/5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% mainEP(flagRun,dataEP,scenario),toc


toc
