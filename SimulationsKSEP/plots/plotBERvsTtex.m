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
 snr=3;
set(0, 'DefaultLineLineWidth', 1.1);


 load BERvsTS1Lim3
 turboAxis=0:dataEP.numberTurbo;
% semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'v:b'),hold on

 SNRdBSoft(snr)
 
 semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'v:b'),hold on
 
 semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'+:','color',mycolor1)
 
 %semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP(:,snr).'],'s:','color',mycolor3),hold on
 semilogy(turboAxis,[BER_soft_PKSEP(snr),BER_soft_PKSEP_turbo(:,snr).'],'s:','color',mycolor2)

 
   load BERvsTS3Lim3
   turboAxis=0:dataEP.numberTurbo;
% semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o--b'),hold on
 semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'+-.','color',mycolor1),
 semilogy(turboAxis,[BER_soft_PKSEP(snr),BER_soft_PKSEP_turbo(:,snr).'],'s-.','color',mycolor2),
 %semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s--c')



%  
%   load BERvsTS3Lim3
%    turboAxis=0:dataEP.numberTurbo;
%  %semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o--b'),hold on
% % semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'s--','color',mycolor3),
%  semilogy(turboAxis,[BER_soft_PKSEP(snr),BER_soft_PKSEP_turbo(:,snr).'],'o-.','color',mycolor2),
%  %semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s--c')
% 
%    load BERvsTS4Lim3
%    turboAxis=0:dataEP.numberTurbo;
%  %semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o--b'),hold on
% % semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'s--','color',mycolor3),
%  semilogy(turboAxis,[BER_soft_PKSEP(snr),BER_soft_PKSEP_turbo(:,snr).'],'d--','color',mycolor2),
%  %semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s--c')
%  
%    load BERvsTS10Lim3
%  turboAxis=0:dataEP.numberTurbo;
% % semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o-b'),hold on
% % semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'s-','color',mycolor3),
%  semilogy(turboAxis,[BER_soft_PKSEP(snr),BER_soft_PKSEP_turbo(:,snr).'],'x-','color',mycolor2),
%  % semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s-.c')
 
 
 legend('KS','KSEP-S1',...
      'KSEP-S3','Location','northeast')
%  legend('MMSE-S1','P-BEP-S1','D-BEP-S1','Sun-S1','MMSE-S10','P-BEP-S10','D-BEP-S10','Sun-S10')
%  title([num2str(dataEP.numberChannels),' channels with ',num2str(dataEP.numberTaps),...
%                 ' taps, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
%                 ', ',
%                 num2str(dataEP.numberFrames*dataEP.numberSimulations),...S
%                 ' frames, ', constelation_name, ' symbols and SNRdB=', num2str(SNRdBSoft(snr)), 'dB']);
%title([ 'limLLR=3 and SNRdB=', num2str(SNRdBSoft(snr)), 'dB']);

grid on
xlabel('T')%('E_b/N_0')
ylabel('BER')
%axis([0 10 1e-4 .5])
set(gca,'FontName','Times','FontSize',14)
%matlab2tikz('BERvsEbN64QAMPRoakisKSEPrngG4R.tex','height','\figureheight','width','\figurewidth')
print('BERvsEbN64QAMPRoakisKSEPrngG4R1','-dpdf')