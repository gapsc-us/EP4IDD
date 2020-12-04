 path(path,'matlab2tikz/src')
 
 
% \definecolor{mycolor1}{rgb}{1.00000,0.50000,0.00000}%
% \definecolor{mycolor2}{rgb}{0.00000,0.80000,1.00000}%
% \definecolor{mycolor3}{rgb}{1.00000,0.00000,1.00000}%
% \definecolor{mycolor4}{rgb}{0.45, 0.31, 0.59}%
% \definecolor{mycolor5}{rgb}{0.6, 0.4, 0.8}



 clear all, %close all
 mycolor=[.1,.5,.1];
 mycolor1=[1,.5,0];
 mycolor2=[0,.8,1];
mycolor3=[1,0,1];
  figure,
 %idx=1, snr=23
 %snr=1;
 kk=5;
set(0, 'DefaultLineLineWidth', 1.1);
load ../Results/ResultsN/BERN256Lim3
 
 %semilogy(SNRdBSoft-3,BER_hard_MMSE,'v-k'),hold on
 set(0, 'DefaultLineLineWidth', 0.5);
 semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v-b'),hold on
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d-','color',mycolor1)
 semilogy(SNRdBSoft,BER_soft_SEP_turbo(kk,:),'s-','color',mycolor)
 semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o-','color',mycolor2)
 set(0, 'DefaultLineLineWidth', 1.1);
 
 load ../Results/ResultsN/BERN512Lim3
  semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v-.b')
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d-.','color',mycolor1)
  semilogy(SNRdBSoft,BER_soft_SEP_turbo(kk,:),'s-.','color',mycolor)
 semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o-.','color',mycolor2)

  load ../Results/ResultsN/BERN1024Lim3
  semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v--b')
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d--','color',mycolor1)
  semilogy(SNRdBSoft,BER_soft_SEP_turbo(kk,:),'s--','color',mycolor)
 semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o--','color',mycolor2)
 
  load ../Results/ResultsN/BERN2048Lim3
  semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v:b')
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d:','color',mycolor1)
  semilogy(SNRdBSoft,BER_soft_SEP_turbo(kk,:),'s:','color',mycolor)
 semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o:','color',mycolor2)

 
   load ../Results/ResultsN/BERN4096Lim3
   set(0, 'DefaultLineLineWidth', 1.5);
  semilogy(SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v-b')
 semilogy(SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'d-','color',mycolor1)
  semilogy(SNRdBSoft,BER_soft_SEP_turbo(kk,:),'s-','color',mycolor)
 semilogy(SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'o-','color',mycolor2)
 
 
%   legend('LMMSE N256','BPEP N256', 'SEP N256', 'KSEP N256',...
%             'LMMSE N512','BPEP N512', 'SEP N512', 'KSEP N512',...
%     'LMMSE N1024','BPEP N1024', 'SEP N1024', 'KSEP N1024',...
%   'LMMSE N2048','BPEP N2048', 'SEP N2048', 'KSEP N2048',...
%       'LMMSE N4096','BPEP N4096', 'SEP N4096', 'KSEP N4096',...
%       'Location','Southwest')%, ...

 
  legend('LMMSE','BP-EP', 'SEP', 'KSEP',...
      'Location','Southwest')%, ...
  

grid on
xlabel('E_b/N_0')
ylabel('BER')

if dataEP.scenario==1,
    pbm='Eq';
    folderName='../figuresEq/';
elseif dataEP.scenario==2,
    pbm='MIMO';
    folderName='../figuresMIMO/';
end


%axis([26 42 3e-5 .5])
set(gca,'FontName','Times','FontSize',14)
if dataEP.flagPSK, modName='PSK'; else modName='QAM'; end
nameFig=[folderName,pbm,'BERvsN',num2str(dataEP.M),modName,dataEP.channelName,'Fr',num2str(dataEP.numberFrames*dataEP.numberSimulations)];
%matlab2tikz(nameFig,'height','\figureheight','width','\figurewidth')
matlab2tikz([nameFig,'.tex'],'height','\figureheight','width','\figurewidth')


title([num2str(dataEP.M),modName,' ', dataEP.channelName])
print(nameFig,'-dpdf')