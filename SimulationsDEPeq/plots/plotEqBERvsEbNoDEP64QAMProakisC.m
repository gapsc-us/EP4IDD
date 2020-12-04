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

set(0, 'DefaultLineLineWidth', 1.1);
load ../ResultsEQ/Results64QAMLim3ProakisC/BERvsSNRdB64QAMLim3ProakisCFr10000

 kk=10;
 semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v-b'),hold on
 
   kk=3;
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'+:r','color',mycolor1)
 
   kk=10;
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'+-r','color',mycolor1)
 
 
 kk=3;
 semilogy(SNRdBSoft,BER_soft_PBEP_turbo(kk,:),'s:m')%,'color',mycolor)
 
 kk=10;
 semilogy(SNRdBSoft,BER_soft_PBEP_turbo(kk,:),'s-m')%,'color',mycolor)
 
%   kk=10;
%  semilogy(SNRdBSoft,BER_soft_PBEP_turbo(kk,:),'+--m')
 

 
%  kk=10;
%  semilogy(SNRdBSoft,BER_soft_DBEP_turbo(kk,:),'+-c')
%   kk=5;
%  semilogy(SNRdBSoft,BER_soft_PFEP_turbo(kk,:),'d:r')
 kk=10;
 semilogy(SNRdBSoft,BER_soft_PFEP_turbo(kk,:),'d-r')

  kk=10;
 semilogy(SNRdBSoft,BER_soft_DFEP_turbo(kk,:),'x-','color',mycolor)
 
 
  kk=3;
 semilogy(SNRdBSoft,BER_soft_DBEP_turbo(kk,:),'o:','color',mycolor2)
 
 kk=10;
 semilogy(SNRdBSoft,BER_soft_DBEP_turbo(kk,:),'o-','color',mycolor2)
 
%     kk=10;
%  semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'+:r')
 
%T=3,5,10 legend('MMSE', 'PBEP-3','PBEP-5','PBEP-10','DBEP-3','DBEP-5','DBEP-10','PFEP-3','BPEP-3','BPEP-5','BPEP-10')

 legend('LMMSE-10','BP-EP-3','BP-EP-10','BEP-3','BEP-10','FEP-10','D-FEP-10','D-BEP-3','D-BEP-10','Location','southwest')













grid on
%axis([8 20 .5e-4 .2])
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