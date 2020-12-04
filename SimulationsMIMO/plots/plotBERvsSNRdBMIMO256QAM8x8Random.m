path(path,'matlab2tikz/src')
 
 
% \definecolor{mycolor1}{rgb}{1.00000,0.50000,0.00000}%
% \definecolor{mycolor2}{rgb}{0.00000,0.80000,1.00000}%
% \definecolor{mycolor3}{rgb}{1.00000,0.00000,1.00000}%
% \definecolor{mycolor4}{rgb}{0.45, 0.31, 0.59}%
% \definecolor{mycolor5}{rgb}{0.6, 0.4, 0.8}



 clear all,% close all
 mycolor=[.1,.5,.1];
 mycolor1=[1,.5,0];
 mycolor2=[0,.8,1];
 mycolor3=[1,0,1];
  figure,
 %idx=1, snr=23
 %snr=1;

set(0, 'DefaultLineLineWidth', 1.1);
%load ../ResultsMIMO/Results256QAMLim3Random/BERvsSNRdB256QAMMIMO128x128Lim3RandomFr100
load ../ResultsMIMO/Results256QAMLim3Random/BERvsSNRdB256QAMMIMO8x8Lim3RandomFr100


%                'SNRdBSoft','dataEP',...


 kk=5;
 semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v-b'),hold on
 % semilogy(SNRdBSoft,BER_soft_MMSEAproxGS_turbo(kk,:),'^-.b') %,'color',mycolor)
  semilogy(SNRdBSoft,BER_soft_EPICLMMSE_turbo(kk,:),'d-','color',mycolor1)
   semilogy(SNRdBSoft,BER_soft_PBEP_turbo(kk,:),'x-','color',mycolor)
  semilogy(SNRdBSoft,BER_soft_DBEP_turbo(kk,:),'s-','color',mycolor2)

% semilogy(SNRdBSoft,BER_soft_DBEPAproxGS_turbo(kk,:),'o-.','color',mycolor2)

% semilogy(SNRdBSoft,BER_soft_DBEPAproxPCG_turbo(kk,:),'+--','color',mycolor2)
% semilogy(SNRdBSoft,BER_soft_DBEPAproxNeuman_turbo(kk,:),'p:','color',mycolor2), 
 
 
  legend('LMMSE','EP-IC','BEP','DEP',...%'DEP-GS',...
  'Location','SouthWest')%, ...
  %'DEP-CG','DEP-NSE','Location','SouthWest')%, ...

xlabel('E_b/N_0')
ylabel('BER')
%axis([12 15 1e-4 .4])
 
grid on
set(gca,'FontName','Times','FontSize',14)
if dataEP.scenario==1,
    pbm='Eq';
    folderName='../figuresEq/';
elseif dataEP.scenario==2,
    pbm='MIMO';
    folderName='../figuresMIMO/';
end

set(gca,'FontName','Times','FontSize',14)
if dataEP.flagPSK, modName='PSK'; else modName='QAM'; end
if dataEP.scenario==2, pbm=[pbm, num2str(dataEP.numberAntennas(1)),'x',num2str(dataEP.numberAntennas(2))],  end
nameFig=[folderName,pbm,'BERvsEbN0',num2str(dataEP.M),modName,'Lim',num2str(dataEP.LLRlim),dataEP.channelName,'Fr',num2str(dataEP.numberFrames*dataEP.numberSimulations/1000),'k'];
matlab2tikz([nameFig,'.tex'],'height','\figureheight','width','\figurewidth')

title([num2str(dataEP.M),modName,' ', dataEP.channelName])
print([nameFig,'.pdf'],'-dpdf')