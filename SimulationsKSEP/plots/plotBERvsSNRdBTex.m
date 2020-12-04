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
 %snr=1;

set(0, 'DefaultLineLineWidth', 1.1);
load results/ResultsSLim3/BERvsSNRdBLim364QAMRandomFr10k
 
 semilogy(SNRdBSoft+3,BER_hard_MMSE,'v-k'),hold on
 kk=5;
 semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v-b'),hold on
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d-','color',mycolor1)
 semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o-','color',mycolor2)
 
 %k=0 NO TURBO
 semilogy(SNRdBSoft,BER_soft_MMSE,'v:b'),
 semilogy(SNRdBSoft,BER_soft_BPEP,'d:','color',mycolor1)
 semilogy(SNRdBSoft,BER_soft_PKSEP,'o:','color',mycolor2)
  
  kk=3;
  semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v--b')
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d--','color',mycolor1)
 semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o--','color',mycolor2)


  legend('LMMSE T5','BP-EP T5','DKSEP-T5',...
      'LMMSE','BP-EP','DKSEP',...
      'LMMSE T3','BP-EP T3','DKSEP-T3',...
      'Location','Southeast')%, ...
grid on
xlabel('E_b/N_0')
ylabel('BER')
%axis([26 42 3e-5 .5])
set(gca,'FontName','Times','FontSize',14)
if dataEP.flagPSK, modName='PSK'; else modName='QAM'; end
nameFig=['BERvsEbN0',num2str(dataEP.M),modName,dataEP.channelName,'Fr',num2str(dataEP.numberFrames*dataEP.numberSimulations/10),'k'];
%matlab2tikz(nameFig,'height','\figureheight','width','\figurewidth')


title([num2str(dataEP.M),modName,' ', dataEP.channelName])
print(nameFig,'-dpdf')