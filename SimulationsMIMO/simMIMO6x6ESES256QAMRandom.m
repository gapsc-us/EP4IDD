close all,clear all
addpath('../')
addpath('../algorithmsCode')

%Used to generate and compare curves with different LLRlim and EPiterations
tic

flagRun=1;
%run configFiles/EP4S0config64QAMLim5Random
run configFiles/MimoS0config256QAM6x6Lim3Random
disp('S=0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
dataEP.EPiterations=0;
mainEP(flagRun,dataEP,scenario),toc


%run configFiles/EP4Sconfig64QAMLim5Random
run configFiles/MimoSconfig256QAM6x6Lim3Random
disp('S=1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
dataEP.EPiterations=1;
mainEP(flagRun,dataEP,scenario),toc

% disp('S=2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% dataEP.EPiterations=2;
% mainEP(flagRun,dataEP,scenario),toc

disp('S=3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
dataEP.EPiterations=3;
mainEP(flagRun,dataEP,scenario),toc

% disp('S=4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% dataEP.EPiterations=4;
% mainEP(flagRun,dataEP,scenario),toc

disp('S=10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
dataEP.EPiterations=10;
mainEP(flagRun,dataEP,scenario),toc

if flagRun==0
    plotBERvsT
end
toc
