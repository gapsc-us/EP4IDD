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
load ../ResultsEQ/Results64QAMLim3Random/BERvsSNRdB64QAMLim3RandomFr100

 kk=10;
 semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v-b'),hold on
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d-','color',mycolor1)
 semilogy(SNRdBSoft,BER_soft_SEP_turbo(kk,:),'s-','color',mycolor)
 semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o-','color',mycolor2)
 
 %best k
 semilogy(SNRdBSoft,min(BER_soft_BPEP_turbo),'d-.','color',mycolor1)
 
 %k=0
 semilogy(SNRdBSoft,BER_soft_MMSE,'v:b'),
 semilogy(SNRdBSoft,BER_soft_BPEP,'d:','color',mycolor1)
 semilogy(SNRdBSoft,BER_soft_SEP,'s:g','color',mycolor)
 semilogy(SNRdBSoft,BER_soft_PKSEP,'o:','color',mycolor2)
 
  kk=3;
  semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v--b'),
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d--','color',mycolor1)
 semilogy(SNRdBSoft,BER_soft_SEP_turbo(kk,:),'s--','color',mycolor)
 semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o-.','color',mycolor2)

 

  legend('KS','BP-EP','SEP','KSEP','Location','Southwest')%, ...
  
 %    'SEP-0','SEP-5','SEP-10','KSEP-0','KSEP-5','KSEP-10')%,'Location','southwest')
 %legend('MMSE-0','MMSE-10','BPEP-0','BPEP-5','BPEP-10', ...
 %    'SEP-0','SEP-5','SEP-10','KSEP-0','KSEP-5','KSEP-10')%,'Location','southwest')
 
%  turboAxis=0:dataEP.numberTurbo;
%  semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o-b'),hold on
%  semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'v-r'),
%  semilogy(turboAxis,[BER_soft_DBEP(snr),BER_soft_DBEP_turbo(:,snr).'],'d-m')
%  semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s-c')
% %  
% %   load BERvsT4S3snr
% %    turboAxis=0:dataEP.numberTurbo;
% %  semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o--b'),hold on
% %  semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'v--r'),
% %  semilogy(turboAxis,[BER_soft_DBEP(snr),BER_soft_DBEP_turbo(:,snr).'],'d--m'),
% %  semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s--c')
% 
%  
%    load BERvsTS10Lim3
%  turboAxis=0:dataEP.numberTurbo;
%  semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o-.b'),hold on
%  semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'v-.r'),
%  semilogy(turboAxis,[BER_soft_DBEP(snr),BER_soft_DBEP_turbo(:,snr).'],'d-.m'),
%   semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s-.c')
%  
% % legend('MMSE-S1','P-BEP-S1','D-BEP-S1','MMSE-S3','P-BEP-S3','D-BEP-S3','MMSE-S10','P-BEP-S10','D-BEP-S10')
%   legend('MMSE-S1','P-BEP-S1','D-BEP-S1','Sun-S1','MMSE-S10','P-BEP-S10','D-BEP-S10','Sun-S10')
% %  title([num2str(dataEP.numberChannels),' channels with ',num2str(dataEP.numberTaps),...
% %                 ' taps, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
% %                 ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
% %                 ' frames, ', constelation_name, ' symbols and SNRdB=', num2str(SNRdBSoft(snr)), 'dB']);
% title([ 'limLLR=3 and SNRdB=', num2str(SNRdBSoft(snr)), 'dB']);
grid on
axis([8 20 .5e-4 .2])
xlabel('E_b/N_0')
ylabel('BER')

% set(gca,'FontName','Times','FontSize',14)
% matlab2tikz('BERvsEbN064QAMRandomT310.tex','height','\figureheight','width','\figurewidth')
% print('BERvsEbN064QAMRandomT310','-dpdf')


if dataEP.scenario==1,
    pbm='Eq';
    folderName='../figuresEq/';
elseif dataEP.scenario==2,
    pbm='MIMO';
    folderName='../figuresMIMO/';
end

set(gca,'FontName','Times','FontSize',14)
if dataEP.flagPSK, modName='PSK'; else modName='QAM'; end
nameFig=[folderName,pbm,'BERvsEbN0',num2str(dataEP.M),modName,'Lim',num2str(dataEP.LLRlim),dataEP.channelName,'Fr',num2str(dataEP.numberFrames*dataEP.numberSimulations/1000),'k'];
matlab2tikz([nameFig,'.tex'],'height','\figureheight','width','\figurewidth')

title([num2str(dataEP.M),modName,' ', dataEP.channelName])
print([nameFig,'.pdf'],'-dpdf')