 clear all, 
  figure,
 %idx=1, snr=23
 snr=1;

 
 load BERvsTS0Lim3
 turboAxis=0:dataEP.numberTurbo;
 semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o-b'),hold on
 semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'v-r'),
 semilogy(turboAxis,[BER_soft_DBEP(snr),BER_soft_DBEP_turbo(:,snr).'],'d-m')
 semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s-c')
%  
%   load BERvsT4S3snr
%    turboAxis=0:dataEP.numberTurbo;
%  semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o--b'),hold on
%  semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'v--r'),
%  semilogy(turboAxis,[BER_soft_DBEP(snr),BER_soft_DBEP_turbo(:,snr).'],'d--m'),
%  semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s--c')

 
   load BERvsTS1Lim3
 turboAxis=0:dataEP.numberTurbo;
 semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'o-.b'),hold on
 semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'v-.r'),
 semilogy(turboAxis,[BER_soft_DBEP(snr),BER_soft_DBEP_turbo(:,snr).'],'d-.m'),
  semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'s-.c')
 
% legend('MMSE-S1','P-BEP-S1','D-BEP-S1','MMSE-S3','P-BEP-S3','D-BEP-S3','MMSE-S10','P-BEP-S10','D-BEP-S10')
  legend('MMSE-S1','P-BEP-S1','D-BEP-S1','Sun-S1','MMSE-S10','P-BEP-S10','D-BEP-S10','Sun-S10')
%  title([num2str(dataEP.numberChannels),' channels with ',num2str(dataEP.numberTaps),...
%                 ' taps, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
%                 ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
%                 ' frames, ', constelation_name, ' symbols and SNRdB=', num2str(SNRdBSoft(snr)), 'dB']);
title([ 'limLLR=3 and SNRdB=', num2str(SNRdBSoft(snr)), 'dB']);

