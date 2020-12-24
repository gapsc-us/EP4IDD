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
 snr1=3;

set(0, 'DefaultLineLineWidth', 1.1);

load ../resultsEQ/Results64QAMLim3ProakisC/BERvsSNRdB64QAMLim3ProakisCFr10000
  SNRdBSoft(snr1)
 turboAxis=0:dataEP.numberTurbo;
 semilogy(turboAxis,[BER_soft_MMSE(snr1),BER_soft_MMSE_turbo(:,snr1).'],'v-b'),hold on 
 semilogy(turboAxis,[BER_soft_BPEP(snr1),BER_soft_BPEP_turbo(:,snr1).'],'d-','color',mycolor1), 
 semilogy(turboAxis,[BER_soft_SEP(snr1),BER_soft_SEP_turbo(:,snr1).'],'s-','color',mycolor),
 semilogy(turboAxis,[BER_soft_PKSEP(snr1),BER_soft_PKSEP_turbo(:,snr1).'],'o-','color',mycolor2), 


 snr2=10;
 SNRdBSoft(snr2)

 turboAxis=0:dataEP.numberTurbo;
 semilogy(turboAxis,[BER_soft_MMSE(snr2),BER_soft_MMSE_turbo(:,snr2).'],'v--b')
  semilogy(turboAxis,[BER_soft_BPEP(snr2),BER_soft_BPEP_turbo(:,snr2).'],'d--','color',mycolor1), 
 semilogy(turboAxis,[BER_soft_SEP(snr2),BER_soft_SEP_turbo(:,snr2).'],'s--','color',mycolor)
 semilogy(turboAxis,[BER_soft_PKSEP(snr2),BER_soft_PKSEP_turbo(:,snr2).'],'o--','color',mycolor2), hold off

 
 
 legend('KS','BPEP','SEP','KSEP','Location','northeast')
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
axis([0 10 1e-3 0.3])

if dataEP.scenario==1,
    pbm='Eq';
    folderName='../figuresEq/';
elseif dataEP.scenario==2,
    pbm='MIMO';
    folderName='../figuresMIMO/';
end


set(gca,'FontName','Times','FontSize',14)
if dataEP.flagPSK, modName='PSK'; else modName='QAM'; end
nameFig=[folderName,pbm,'BERvsT',num2str(dataEP.M),modName,dataEP.channelName,'Fr',num2str(dataEP.numberFrames*dataEP.numberSimulations/1000),'k'];
matlab2tikz([nameFig,'.tex'],'height','\figureheight','width','\figurewidth')

title([num2str(dataEP.M),modName,' ', dataEP.channelName,' for SNR=',num2str(SNRdBSoft(snr1)),', ',num2str(SNRdBSoft(snr2))])
print([nameFig,'.pdf'],'-dpdf')

%axis([0 10 1e-4 .5])
% set(gca,'FontName','Times','FontSize',14)
% matlab2tikz('BERvsEbN64QAMPRoakisKSEPrngG4R.tex','height','\figureheight','width','\figurewidth')
% print('BERvsEbN64QAMPRoakisKSEPrngG4R1','-dpdf')

