path(path,'matlab2tikz/src')


% \definecolor{mycolor1}{rgb}{1.00000,0.50000,0.00000}%
% \definecolor{mycolor2}{rgb}{0.00000,0.80000,1.00000}%
% \definecolor{mycolor3}{rgb}{1.00000,0.00000,1.00000}%
% \definecolor{mycolor4}{rgb}{0.45, 0.31, 0.59}%
% \definecolor{mycolor5}{rgb}{0.6, 0.4, 0.8}



clear all, close all
mycolor=[.1,.5,.1];
mycolor1=[1,.5,0];
mycolor2=[0,.8,1];
mycolor3=[1,0,1];
figure,
%idx=1, snr=23
snr=2;


%  load BERvsTS0Lim3

set(0, 'DefaultLineLineWidth', 1.1);
%load ../ResultsMIMO/Results256QAMLim3Random/BERvsSNRdB256QAMMIMO128x128Lim3RandomFr100
load ../ResultsMIMO/Results256QAMLim3RandomS0/BERvsTS0256QAMMIMO6x6Lim3RandomFr100



turboAxis=0:dataEP.numberTurbo;
% semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o:b')
semilogy(turboAxis,[BER_soft_EPICLMMSE(snr),BER_soft_EPICLMMSE_turbo(:,snr).'],'d-','color',mycolor1),hold on
semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'x-','color',mycolor),
semilogy(turboAxis,[BER_soft_DBEP(snr),BER_soft_DBEP_turbo(:,snr).'],'s-','color',mycolor2),
% semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s-c')


load ../ResultsMIMO/Results256QAMLim3RandomS1/BERvsTS1256QAMMIMO6x6Lim3RandomFr100
turboAxis=0:dataEP.numberTurbo;
% semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o-b'),hold on
semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'x-.','color',mycolor),
semilogy(turboAxis,[BER_soft_DBEP(snr),BER_soft_DBEP_turbo(:,snr).'],'s-.','color',mycolor2)
% semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s-c')
%


load ../ResultsMIMO/Results256QAMLim3RandomS3/BERvsTS3256QAMMIMO6x6Lim3RandomFr100
turboAxis=0:dataEP.numberTurbo;
% semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o--b'),hold on
semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'x--','color',mycolor),
semilogy(turboAxis,[BER_soft_DBEP(snr),BER_soft_DBEP_turbo(:,snr).'],'s--','color',mycolor2),
%semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s--c')


load ../ResultsMIMO/Results256QAMLim3RandomS10/BERvsTS10256QAMMIMO6x6Lim3RandomFr100
turboAxis=0:dataEP.numberTurbo;
% semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o-.b'),hold on
semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'x:','color',mycolor),
semilogy(turboAxis,[BER_soft_DBEP(snr),BER_soft_DBEP_turbo(:,snr).'],'s:','color',mycolor2),


%NOTA: 'LMMSE', %=p-BEP S0
legend('EP-IC','BEP S0','DEP S0','BEP S1','DEP S1',...
    'BEP S3','DEP S3','BEP S10','DEP S10','Location','southwest')


grid on
xlabel('T')%('E_b/N_0')
ylabel('BER')
%axis([0 10 3e-4 .05])
set(gca,'FontName','Times','FontSize',14)
matlab2tikz('BERvsT256QAM6x6T10.tex','height','\figureheight','width','\figurewidth')
title(['256QAM, 6x6,T=10, E_b/N_0=',num2str(SNRdBSoft(snr))])
print('BERvsT256QAM6x6T10','-dpdf')
