close all,clear all
addpath('../')
addpath('../algorithmsCode')
%addpath('../ParityCheckMatrix')
%Used to generate and compare curves with different LLRlim and EPiterations
tic
flagRun=1;

run configFiles/EP4Nconfig16QAMLim3ProakisC
disp('N=256 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
dataEP.channelBlockLength=256;dataEP.bits=log2(dataEP.channelBlockLength)-1;
mainEP(flagRun,dataEP,scenario),toc

disp('N=512 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
dataEP.channelBlockLength=512;dataEP.bits=log2(dataEP.channelBlockLength)-1;
mainEP(flagRun,dataEP,scenario),toc

disp('N=1024 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
dataEP.channelBlockLength=1024;dataEP.bits=log2(dataEP.channelBlockLength)-1;
mainEP(flagRun,dataEP,scenario),toc

disp('N=2048 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
dataEP.channelBlockLength=2048;dataEP.bits=log2(dataEP.channelBlockLength)-1;
mainEP(flagRun,dataEP,scenario),toc

disp('N=4096 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
dataEP.channelBlockLength=4096;dataEP.bits=log2(dataEP.channelBlockLength)-1;
mainEP(flagRun,dataEP,scenario),toc
% 
% disp('S=2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% data.EPiterations=2;
% mainEP(flagRun,data),toc
% 
% disp('S=3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% data.EPiterations=3;
% mainEP(flagRun,data),toc
% 
% disp('S=4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% data.EPiterations=4;
% mainEP(flagRun,data),toc
% 
% disp('S=10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% data.EPiterations=10;
% mainEP(flagRun,data),toc

if flagRun==0
    plotBERvsT
end
toc
