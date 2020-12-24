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
load ../ResultsEQ/Results256QAMLim3Random/BERvsSNRdB256QAMLim3RandomFr100.mat
kk=10;
 semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v-b'),hold on
 kk=10;
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d-','color',mycolor1)
 kk=10;
 semilogy(SNRdBSoft,BER_soft_SEP_turbo(kk,:),'s-','color',mycolor)
  kk=10;
 semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o-','color',mycolor2)
 
semilogy(SNRdBSoft,min(BER_soft_BPEP_turbo),'d-.','color',mycolor1)
 
 semilogy(SNRdBSoft,BER_soft_MMSE,'v:b'),hold on
%   kk=5;
%  semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v-.b'),hold on


  
 %kk=3;
 %semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d:','color',mycolor1)
 semilogy(SNRdBSoft,BER_soft_BPEP,'d:','color',mycolor1)
  kk=3;
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d--','color',mycolor1)

 
%  kk=3;
%  semilogy(SNRdBSoft,BER_soft_SEP_turbo(kk,:),'s:g','color',mycolor)
 semilogy(SNRdBSoft,BER_soft_SEP,'s:g','color',mycolor)
 kk=3;
 semilogy(SNRdBSoft,BER_soft_SEP_turbo(kk,:),'s--','color',mycolor)

 
 %kk=3;
 %semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o:','color',mycolor2)
 semilogy(SNRdBSoft,BER_soft_PKSEP,'o:','color',mycolor2)
 kk=3;
 semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o-.','color',mycolor2)

 

  legend('KS','BP-EP','SEP','KSEP','Location','Southwest')%, ...
grid on
xlabel('E_b/N_0')
ylabel('BER')
axis([12 20 1e-5 .3])
% set(gca,'FontName','Times','FontSize',14)
% matlab2tikz('BERvsEbN256QAMProakisT310.tex','height','\figureheight','width','\figurewidth')
% title('256QAM, Proakis C')
% print('BERvsEbN0256QAMProakisCT310','-dpdf')


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