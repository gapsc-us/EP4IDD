function plotEPeqSNRdB(dataEP,saveName,scenario)
%
% plotEPeqSNRdB(dataEP,saveName,scenario)
%
% Author: Irene Santos Vel?zquez
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 23/05/2018
%
% Description: This function plots the results that have been written in a
% log file. The simulation stores values for several runs in parallel. This
% function reads them all and generates a .mat file with the values of the curves
% (averaged over the different runs)
%
% Inputs:
%  dataEP is a structure with several fields such as:
%     dataEP.numberFrames is the number of words transmitted per simulation
%     n is the length of the transmitted block in symbols.
%     dataEP.bits is the size of information dataEP.bits
%     dataEP.numberTaps is the number of taps of the channel in equalization
%     dataEP.numberAntennas is the number of transmitted and receive antennas in MIMO
%     dataEP.numberChannels is the number of simulated channels
%     dataEP.M is the constellation size
%     dataEP.flagPSK is 1 if the complex constellation is a PSK and 0 if it is a QAM
%     dataEP.flagIterDec is 1 if the figures of number of iterations of the channel
%       decoder is plotted
%  saveName indicates the name of file to save the .mat with the overal results
%  scenario if 1 equalization is executed, if 2 MIMO is executed

DELIMITER=' ';
HEADERLINES=4;

n=dataEP.channelBlockLength;

if ~exist('dataEP.folderSave'), dataEp.folderSave='Results'; end
[status, msg, msgID] = mkdir(dataEP.folderSave);

if ~exist('dataEP.rate'), rate=1/2 ;%Default rate
else
rate=dataEP.rate; 
end

if dataEP.scenario==1 % equalization
    name=([dataEP.folderSave,'/','resultsEPeq', '_b' num2str(dataEP.bits) '_M' num2str(dataEP.M) '_L' num2str(dataEP.numberTaps) '_n' num2str(n) '_numCh' num2str(dataEP.numberChannels) '_' ]);
else
    name=([dataEP.folderSave,'/','resultsEPmm', '_b' num2str(dataEP.bits) '_M' num2str(dataEP.M) '_NtxNr' num2str(dataEP.numberAntennas(1)) 'x' num2str(dataEP.numberAntennas(2)) '_n' num2str(n) '_numCh' num2str(dataEP.numberChannels) '_' ]);
end

% Checking the SNRdBSim, dataEP.numberTurbo, dataEP.numberTaps and dataEP.numberChannels (because
% in the main file we ignore theirs values if the are changed inside
% blockEPSNRdB_cluster function
nameOmain=([name num2str(1) '.txt']);
fidauxMain=fopen(nameOmain,'r');
newData1=importdata(nameOmain,DELIMITER, HEADERLINES);
TextData=newData1.textdata;
Data=str2num(TextData{3});
SNRdB=Data(1):Data(2):Data(3);
dataEP.numberTurbo=Data(4);
dataEP.numberChannels=Data(5);
if dataEP.scenario==1
    dataEP.numberTaps=Data(6);
else
    dataEP.numberAntennas=[Data(6) Data(7)];
end
fclose(fidauxMain);

SNRdBSoft=SNRdB+3;
if dataEP.flagIterDec
    numberMethods=53+53*dataEP.numberTurbo;
else
    numberMethods=34+34*dataEP.numberTurbo;
end
numberPointsSNRdB=length(SNRdB);

berChannelsSimulations=zeros(1,numberPointsSNRdB*numberMethods);
dataEP.numberSimulations, dataEP.numberChannels
for simulationNumber=1:dataEP.numberSimulations
    nameOmain=([name num2str(simulationNumber) '.txt'])
    
    fidauxMain=fopen(nameOmain,'r');
    %Read the data from the file and estimate the errores
    newData1=importdata(nameOmain,DELIMITER, HEADERLINES);
    
    berMatrix=newData1.data;
    berMatrix=berMatrix(:,2:end); %Firsr column is the number of channel
    
    berChannels=sum(berMatrix(1:dataEP.numberChannels,:),1)/dataEP.numberChannels;
    berChannelsSimulations=berChannelsSimulations+berChannels;
    
    fclose(fidauxMain);
end
berChannelsSimulations=berChannelsSimulations/dataEP.numberSimulations;

berChannelsSimulations=reshape(berChannelsSimulations,numberMethods,numberPointsSNRdB);

indnext=1;
% We go writing into the file...
BER_hard_MMSE=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_MMSE_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_MMSEAproxGS=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_MMSEAproxGS_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_MMSEAproxNeuman=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_MMSEAproxNeuman_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_MMSEAproxPCG=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_MMSEAproxPCG_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_FMMSE=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_FMMSE_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_BEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_BEP_turbo(k1-1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_PBEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_PBEP_turbo(k1-1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_DBEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_DBEP_turbo(k1-1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_DBEPAproxGS=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_DBEPAproxGS_turbo(k1-1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_DBEPAproxNeuman=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_DBEPAproxNeuman_turbo(k1-1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_DBEPAproxPCG=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_DBEPAproxPCG_turbo(k1-1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_PFEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_PFEP_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_SEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_SEP_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_EPICLMMSE=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_EPICLMMSE_turbo(k1-1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_hard_Optimal=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_hard_Optimal_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end

BER_soft_MMSE=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_MMSE_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_MMSEAproxGS=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_MMSEAproxGS_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_MMSEAproxNeuman=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_MMSEAproxNeuman_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_MMSEAproxPCG=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_MMSEAproxPCG_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_FMMSE=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_FMMSE_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_BEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_BEP_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_PBEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_PBEP_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_DBEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_DBEP_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_DBEPAproxGS=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_DBEPAproxGS_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_DBEPAproxNeuman=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_DBEPAproxNeuman_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_DBEPAproxPCG=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_DBEPAproxPCG_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_PFEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_PFEP_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_DFEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_DFEP_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_SEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_SEP_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_PKSEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_PKSEP_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_DKSEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_DKSEP_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_BPEP=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_BPEP_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_EPICLMMSE=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_EPICLMMSE_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end
BER_soft_Optimal=berChannelsSimulations(indnext,:);
indini=indnext+1;
indnext=indini;
for k1=indini:indini+dataEP.numberTurbo-1
    BER_soft_Optimal_turbo(k1-indini+1,:)=berChannelsSimulations(k1,:);
    indnext=k1+1;
end


%% IterDec
if dataEP.flagIterDec
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_MMSEAproxGS(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_MMSEAproxNeuman(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_MMSEAproxPCG(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_MMSE(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_FMMSE(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_BEP(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_PBEP(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_DBEP(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_DBEPAproxGS(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_DBEPAproxNeuman(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_DBEPAproxPCG(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_PFEP(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_DFEP(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_SEP(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_PKSEP(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_DKSEP(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_EPICLMMSE(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_BPEP(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
    indini=indnext;
    for k1=indini:indini+dataEP.numberTurbo
        niterdec_Optimal(k1-indini+1,:)=berChannelsSimulations(k1,:);
        indnext=k1+1;
    end
end

if (dataEP.complexFlag==0 && dataEP.M==2)
    constelation_name='BPSK';
elseif (dataEP.complexFlag==0 && dataEP.M>2)
    constelation_name=[num2str(dataEP.M),'-PAM'];
elseif (dataEP.complexFlag==1)
    if dataEP.flagPSK
        constelation_name=[num2str(dataEP.M),'-PSK'];
    else
        constelation_name=[num2str(dataEP.M),'-QAM'];
    end
end

%% NO TURBO
if dataEP.numberTurbo==0
    % Hard
    flagBEP=scenario.flagBEP; %(sum(BER_hard_BEP));
    flagMMSE=scenario.flagMMSE; %(sum(BER_hard_MMSE));
    flagSEP=scenario.flagSEP; %(sum(BER_hard_SEP));
    flagOptimal=scenario.flagOptimal; %(sum(BER_hard_Optimal));
    flagEPICLMMSE=scenario.flagEPICLMMSE; %sum(BER_hard_EPICLMMSE);
    flagPKSEP=scenario.flagPKSEP ;%sum(sum(BER_soft_PKSEP_turbo));
    flagPBEP=scenario.flagPBEP; %sum(BER_hard_PBEP);
    flagDBEP=scenario.flagDBEP; %sum(BER_hard_DBEP);
    flagPFEP=scenario.flagPFEP; %sum(BER_hard_PFEP);
    flagFMMSE=scenario.flagFMMSE; %sum(BER_hard_FMMSE);
    
    %% No turbo & Hard
    if ~dataEP.flagOnlySoft
        figure
        index=0;
        if flagMMSE, semilogy(SNRdBSoft,BER_hard_MMSE,'v-b'),hold on,...
                index=index+1; legendText(index)={['MMSE No-T']};  end
        if flagFMMSE, semilogy(SNRdBSoft,BER_hard_FMMSE,'o-m'),hold on,...
                index=index+1; legendText(index)={['FMMSE No-T']};  end
        if flagBEP, semilogy(SNRdBSoft,BER_hard_BEP,'s-r'),hold on,...
                index=index+1; legendText(index)={['BEP No-T']};  end
        if flagPBEP, semilogy(SNRdBSoft,BER_hard_PBEP,'v-r'),hold on,...
                index=index+1; legendText(index)={['P-BEP No-T']};  end
        if flagDBEP, semilogy(SNRdBSoft,BER_hard_DBEP,'d-r'),hold on,...
                index=index+1; legendText(index)={['D-BEP No-T']};  end
        if flagPFEP, semilogy(SNRdBSoft,BER_hard_PFEP,'d-c'),hold on,...
                index=index+1; legendText(index)={['P-FEP No-T']};  end
        if flagSEP, semilogy(SNRdBSoft,BER_hard_SEP,'o-g'),hold on,...
                index=index+1; legendText(index)={['SEP No-T']};  end
        if flagEPICLMMSE, semilogy(SNRdBSoft,BER_hard_EPICLMMSE,'+-m'),hold on,...
                index=index+1; legendText(index)={['EP IC-LMMSE No-T']};   end
        if flagOptimal, semilogy(SNRdBSoft,BER_hard_Optimal,'d-k'),hold on,...
                index=index+1; legendText(index)={['Optimal No-T']};  end
        
        hold off
        legend(legendText)
        xlabel('SNRdB (dB)'),...
            ylabel('BER hard')
        
        if dataEP.scenario==1
            title([num2str(dataEP.numberChannels),' channels with ',num2str(dataEP.numberTaps),...
                ' taps, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
                ', ', num2str(dataEP.dataEP.numberFrames*dataEP.numberSimulations),...
                ' frames and ', constelation_name, ' symbols']);
        else
            title([num2str(dataEP.numberChannels),' channels with Ntx=',num2str(dataEP.numberAntennas(1)),...
                ' and Nrx=',num2str(dataEP.numberAntennas(2)),...
                ' antennas, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
                ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
                ' frames and ', constelation_name, ' symbols']);
        end
        
        
    end
    
    %% No turbo & Soft
    
    figure
    index=0;
    if flagMMSE, semilogy(SNRdBSoft,BER_soft_MMSE,'v-b'),hold on,...
            index=index+1; legendText(index)={['MMSE No-T']};  end
    if flagFMMSE, semilogy(SNRdBSoft,BER_soft_FMMSE,'o-m'),hold on,...
            index=index+1; legendText(index)={['FMMSE No-T']};  end
    if flagBEP, semilogy(SNRdBSoft,BER_soft_BEP,'s-r'),hold on,...
            index=index+1; legendText(index)={['BEP No-T']};  end
    if flagPBEP, semilogy(SNRdBSoft,BER_soft_PBEP,'v-r'),hold on,...
            index=index+1; legendText(index)={['P-BEP No-T']};  end
    if flagDBEP, semilogy(SNRdBSoft,BER_soft_DBEP,'d-r'),hold on,...
            index=index+1; legendText(index)={['D-BEP No-T']};  end
    if flagPFEP, semilogy(SNRdBSoft,BER_soft_PFEP,'d-c'),hold on,...
            index=index+1; legendText(index)={['P-FEP No-T']};  end
    if flagSEP, semilogy(SNRdBSoft,BER_soft_SEP,'o-g'),hold on,...
            index=index+1; legendText(index)={['SEP No-T']};  end
    if flagPKSEP, semilogy(SNRdBSoft,BER_soft_PKSEP,'*-g'),hold on,...
            index=index+1; legendText(index)={['P-KSEP No-T']}; end
    if flagEPICLMMSE, semilogy(SNRdBSoft,BER_soft_EPICLMMSE,'+-m'),hold on,...
            index=index+1; legendText(index)={['EP IC-LMMSE No-T']};   end
    if flagOptimal, semilogy(SNRdBSoft,BER_soft_Optimal,'d-k'),hold on,...
            index=index+1; legendText(index)={['Optimal No-T']};  end
    
    hold off
    legend(legendText)
    xlabel('SNRdB (dB)'),...
        ylabel('BER soft')
    if dataEP.scenario==1
        title([num2str(dataEP.numberChannels),' channels with ',num2str(dataEP.numberTaps),...
            ' taps, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
            ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
            ' frames and ', constelation_name, ' symbols']);
    else
        title([num2str(dataEP.numberChannels),' channels with Ntx=',num2str(dataEP.numberAntennas(1)),...
            ' and Nrx=',num2str(dataEP.numberAntennas(2)),...
            ' antennas, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
            ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
            ' frames and ', constelation_name, ' symbols']);
    end
    axis([SNRdBSoft(1) SNRdBSoft(end) 1e-6 1])
    
    
    %% TURBO
else
    flagOptimal=sum(sum(BER_hard_Optimal_turbo));
    flagSEP=sum(sum(BER_soft_SEP_turbo));
    flagBEP=sum(sum(BER_hard_BEP_turbo));
    flagEPICLMMSE=sum(sum(BER_hard_EPICLMMSE_turbo));
    if exist('BER_soft_BPEP_turbo'), flagBEPSun=sum(sum(BER_soft_BPEP_turbo)); else flagBEPSun=0; end
    flagPKSEP=sum(sum(BER_soft_PKSEP_turbo));
    flagDKSEP=sum(sum(BER_soft_DKSEP_turbo));
    flagPBEP=sum(sum(BER_hard_PBEP_turbo));
    flagDBEP=sum(sum(BER_hard_DBEP_turbo));
    flagDBEPAproxGS=sum(sum(BER_hard_DBEPAproxGS_turbo));
    flagDBEPAproxNeuman=sum(sum(BER_hard_DBEPAproxNeuman_turbo));
    flagDBEPAproxPCG=sum(sum(BER_hard_DBEPAproxPCG_turbo));
    flagPFEP=sum(sum(BER_hard_PFEP_turbo));
    flagDFEP=sum(sum(BER_soft_DFEP_turbo));
    flagMMSE=sum(sum(BER_hard_MMSE_turbo));
    flagMMSEAproxGS=sum(sum(BER_hard_MMSEAproxGS_turbo));
    flagMMSEAproxNeuman=sum(sum(BER_hard_MMSEAproxNeuman_turbo));
    flagMMSEAproxPCG=sum(sum(BER_hard_MMSEAproxPCG_turbo));
    flagFMMSE=sum(sum(BER_hard_FMMSE_turbo));
    
    if ~dataEP.flagOnlySoft %Turbo & Hard
        
        for kk=1:dataEP.numberTurbo
            % Hard
            figure
            index=0;
            if flagMMSE, semilogy(SNRdBSoft,BER_hard_MMSE,'v-b',SNRdBSoft,BER_hard_MMSE_turbo(kk,:),'v--b'),hold on,...
                    index=index+1; legendText(index)={['MMSE No-T']}; index=index+1; legendText(index)={['MMSE T',num2str(kk)]};  end
            if flagFMMSE, semilogy(SNRdBSoft,BER_hard_FMMSE,'o-m',SNRdBSoft,BER_hard_FMMSE_turbo(kk,:),'o--m'),hold on,...
                    index=index+1; legendText(index)={['FMMSE No-T']}; index=index+1; legendText(index)={['FMMSE T',num2str(kk)]};  end
            if flagBEP, semilogy(SNRdBSoft,BER_hard_BEP,'s-r',SNRdBSoft,BER_hard_BEP_turbo(kk,:),'s--r'),hold on,...
                    index=index+1; legendText(index)={['BEP No-T']}; index=index+1; legendText(index)={['BEP T',num2str(kk)]};  end
            if flagPBEP, semilogy(SNRdBSoft,BER_hard_PBEP,'v-r',SNRdBSoft,BER_hard_PBEP_turbo(kk,:),'v--r'),hold on,...
                    index=index+1; legendText(index)={['P-BEP No-T']}; index=index+1; legendText(index)={['P-BEP T',num2str(kk)]};  end
            if flagDBEP, semilogy(SNRdBSoft,BER_hard_DBEP,'d-r',SNRdBSoft,BER_hard_DBEP_turbo(kk,:),'d--r'),hold on,...
                    index=index+1; legendText(index)={['D-BEP No-T']}; index=index+1; legendText(index)={['D-BEP T',num2str(kk)]};  end
            if flagPFEP, semilogy(SNRdBSoft,BER_hard_PFEP,'d-c',SNRdBSoft,BER_hard_PFEP_turbo(kk,:),'d--c'),hold on,...
                    index=index+1; legendText(index)={['P-FEP No-T']}; index=index+1; legendText(index)={['P-FEP T',num2str(kk)]};  end
            if flagSEP, semilogy(SNRdBSoft,BER_hard_SEP,'o-g',SNRdBSoft,BER_hard_SEP_turbo(kk,:),'o--g'),hold on,...
                    index=index+1; legendText(index)={['SEP No-T']}; index=index+1; legendText(index)={['SEP T',num2str(kk)]};  end
            if flagEPICLMMSE, semilogy(SNRdBSoft,BER_hard_EPICLMMSE,'+-m',SNRdBSoft,BER_hard_EPICLMMSE_turbo(kk,:),'+--m'),hold on,...
                    index=index+1; legendText(index)={['EP IC-LMMSE No-T']}; index=index+1; legendText(index)={['EP IC-LMMSE T',num2str(kk)]};  end
            if flagOptimal, semilogy(SNRdBSoft,BER_hard_Optimal,'d-k',SNRdBSoft,BER_hard_Optimal_turbo(kk,:),'d--k'),hold on,...
                    index=index+1; legendText(index)={['Optimal No-T']}; index=index+1; legendText(index)={['Optimal T',num2str(kk)]};  end
            
            hold off
            legend(legendText)
            xlabel('SNRdB (dB)'),...
                ylabel('BER hard')
            
            if dataEP.scenario==1
                title([num2str(dataEP.numberChannels),' channels with ',num2str(dataEP.numberTaps),...
                    ' taps, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
                    ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
                    ' frames and ', constelation_name, ' symbols']);
            else
                title([num2str(dataEP.numberChannels),' channels with Ntx=',num2str(dataEP.numberAntennas(1)),...
                    ' and Nrx=',num2str(dataEP.numberAntennas(2)),...
                    ' antennas, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
                    ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
                    ' frames and ', constelation_name, ' symbols']);
            end
            axis([SNRdB(1) SNRdB(end) 1e-6 1])
            
            
        end
    end
    
    
    %% SOFT
    
    for kk=1:dataEP.numberTurbo
        
        figure
        index=0;
        if flagMMSE, semilogy(SNRdBSoft,BER_soft_MMSE,'v-b',SNRdBSoft,BER_soft_MMSE_turbo(kk,:),'v--b'),hold on,...
                index=index+1; legendText(index)={['MMSE No-T']}; index=index+1; legendText(index)={['MMSE T',num2str(kk)]};  end
        if flagMMSEAproxGS, semilogy(SNRdBSoft,BER_soft_MMSEAproxGS,'o-g',SNRdBSoft,BER_soft_MMSEAproxGS_turbo(kk,:),'o--g'),hold on,...
                index=index+1; legendText(index)={['MMSEaproxGS No-T']}; index=index+1; legendText(index)={['MMSEaproxGS T',num2str(kk)]};  end
        if flagMMSEAproxNeuman, semilogy(SNRdBSoft,BER_soft_MMSEAproxNeuman,'o-m',SNRdBSoft,BER_soft_MMSEAproxNeuman_turbo(kk,:),'o--m'),hold on,...
                index=index+1; legendText(index)={['MMSEaproxNeuman No-T']}; index=index+1; legendText(index)={['MMSEaproxNeuman T',num2str(kk)]};  end
        if flagMMSEAproxPCG, semilogy(SNRdBSoft,BER_soft_MMSEAproxPCG,'o-r',SNRdBSoft,BER_soft_MMSEAproxPCG_turbo(kk,:),'o--r'),hold on,...
                index=index+1; legendText(index)={['MMSEaproxPCG No-T']}; index=index+1; legendText(index)={['MMSEaproxPCGn T',num2str(kk)]};  end
        if flagFMMSE, semilogy(SNRdBSoft,BER_soft_FMMSE,'o-m',SNRdBSoft,BER_soft_FMMSE_turbo(kk,:),'o--m'),hold on,...
                index=index+1; legendText(index)={['FMMSE No-T']}; index=index+1; legendText(index)={['FMMSE T',num2str(kk)]};  end
        if flagBEP, semilogy(SNRdBSoft,BER_soft_BEP,'s-r',SNRdBSoft,BER_soft_BEP_turbo(kk,:),'s--r'),hold on,...
                index=index+1; legendText(index)={['BEP No-T']}; index=index+1; legendText(index)={['BEP T',num2str(kk)]};  end
        if flagPBEP, semilogy(SNRdBSoft,BER_soft_PBEP,'v-r',SNRdBSoft,BER_soft_PBEP_turbo(kk,:),'v--r'),hold on,...
                index=index+1; legendText(index)={['P-BEP No-T']}; index=index+1; legendText(index)={['P-BEP T',num2str(kk)]};  end
        if flagDBEP, semilogy(SNRdBSoft,BER_soft_DBEP,'d-r',SNRdBSoft,BER_soft_DBEP_turbo(kk,:),'d--r'),hold on,...
                index=index+1; legendText(index)={['D-BEP No-T']}; index=index+1; legendText(index)={['D-BEP T',num2str(kk)]};  end
        if flagDBEPAproxGS, semilogy(SNRdBSoft,BER_soft_DBEPAproxGS,'*-c',SNRdBSoft,BER_soft_DBEPAproxGS_turbo(kk,:),'*--c'),hold on,...
                index=index+1; legendText(index)={['D-BEPaproxGS No-T']}; index=index+1; legendText(index)={['D-BEPaproxGS T',num2str(kk)]};  end
        if flagDBEPAproxNeuman, semilogy(SNRdBSoft,BER_soft_DBEPAproxNeuman,'*-k',SNRdBSoft,BER_soft_DBEPAproxNeuman_turbo(kk,:),'*--k'),hold on,...
                index=index+1; legendText(index)={['D-BEPaproxNeuman No-T']}; index=index+1; legendText(index)={['D-BEPaproxNeuman T',num2str(kk)]};  end
        if flagDBEPAproxPCG, semilogy(SNRdBSoft,BER_soft_DBEPAproxPCG,'*-g',SNRdBSoft,BER_soft_DBEPAproxPCG_turbo(kk,:),'*--g'),hold on,...
                index=index+1; legendText(index)={['D-BEPaproxPCG No-T']}; index=index+1; legendText(index)={['D-BEPaproxPCG T',num2str(kk)]};  end
        if flagPFEP, semilogy(SNRdBSoft,BER_soft_PFEP,'d-c',SNRdBSoft,BER_soft_PFEP_turbo(kk,:),'d--c'),hold on,...
                index=index+1; legendText(index)={['P-FEP No-T']}; index=index+1; legendText(index)={['P-FEP T',num2str(kk)]};  end
        if flagDFEP, semilogy(SNRdBSoft,BER_soft_DFEP,'*-c',SNRdBSoft,BER_soft_DFEP_turbo(kk,:),'*--c'),hold on,...
                index=index+1; legendText(index)={['D-FEP No-T']}; index=index+1; legendText(index)={['D-FEP T',num2str(kk)]};  end
        if flagSEP, semilogy(SNRdBSoft,BER_soft_SEP,'o-g',SNRdBSoft,BER_soft_SEP_turbo(kk,:),'o--g'),hold on,...
                index=index+1; legendText(index)={['SEP No-T']}; index=index+1; legendText(index)={['SEP T',num2str(kk)]};  end
        if flagPKSEP, semilogy(SNRdBSoft,BER_soft_PKSEP,'*-g',SNRdBSoft,BER_soft_PKSEP_turbo(kk,:),'*--g'),hold on,...
                index=index+1; legendText(index)={['P-KSEP No-T']}; index=index+1; legendText(index)={['P-KSEP T',num2str(kk)]};  end
        if flagDKSEP, semilogy(SNRdBSoft,BER_soft_DKSEP,'*-r',SNRdBSoft,BER_soft_DKSEP_turbo(kk,:),'*--r'),hold on,...
                index=index+1; legendText(index)={['D-KSEP No-T']}; index=index+1; legendText(index)={['D-KSEP T',num2str(kk)]};  end
        if flagBEPSun, semilogy(SNRdBSoft,BER_soft_BPEP,'+-c',SNRdBSoft,BER_soft_BPEP_turbo(kk,:),'+--c'),hold on,...
                index=index+1; legendText(index)={['BPEP No-T']}; index=index+1; legendText(index)={['BPEP T',num2str(kk)]};  end
        if flagEPICLMMSE, semilogy(SNRdBSoft,BER_soft_EPICLMMSE,'+-m',SNRdBSoft,BER_soft_EPICLMMSE_turbo(kk,:),'+--m'),hold on,...
                index=index+1; legendText(index)={['EP IC-LMMSE No-T']}; index=index+1; legendText(index)={['EP IC-LMMSE T',num2str(kk)]};  end
        if flagOptimal, semilogy(SNRdBSoft,BER_soft_Optimal,'d-k',SNRdBSoft,BER_soft_Optimal_turbo(kk,:),'d--k'),hold on,...
                index=index+1; legendText(index)={['Optimal No-T']}; index=index+1; legendText(index)={['Optimal T',num2str(kk)]};  end
        
        hold off
        legend(legendText)
        xlabel('SNRdB (dB)'),...
            ylabel('BER soft')
        
        
        if dataEP.scenario==1
            title([num2str(dataEP.numberChannels),' channels with ',num2str(dataEP.numberTaps),...
                ' taps, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
                ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
                ' frames and ', constelation_name, ' symbols']);
        else
            title([num2str(dataEP.numberChannels),' channels with Ntx=',num2str(dataEP.numberAntennas(1)),...
                ' and Nrx=',num2str(dataEP.numberAntennas(2)),...
                ' antennas, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
                ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
                ' frames and ', constelation_name, ' symbols']);
        end
        axis([SNRdBSoft(1) SNRdBSoft(end) 1e-6 1])
    end
    
end

if dataEP.flagBERalongT
    turboAxis=0:dataEP.numberTurbo;
    for snr=1:length(SNRdBSoft)
        figure
        index=0;
        clearvars legendText
        
        if flagMMSE, semilogy(turboAxis,[BER_soft_MMSE(snr),BER_soft_MMSE_turbo(:,snr).'],'v-b'),hold on,...
                index=index+1; legendText(index)={['MMSE']}; end
        if flagMMSEAproxGS, semilogy(turboAxis,[BER_soft_MMSEAproxGS(snr),BER_soft_MMSEAproxGS_turbo(:,snr).'],'o-g'),hold on,...
                index=index+1; legendText(index)={['MMSEaproxGS']}; end
        if flagMMSEAproxNeuman, semilogy(turboAxis,[BER_soft_MMSEAproxNeuman(snr),BER_soft_MMSEAproxNeuman_turbo(:,snr).'],'o-m'),hold on,...
                index=index+1; legendText(index)={['MMSEaproxNeuman']}; end
        if flagMMSEAproxPCG, semilogy(turboAxis,[BER_soft_MMSEAproxPCG(snr),BER_soft_MMSEAproxPCG_turbo(:,snr).'],'o-r'),hold on,...
                index=index+1; legendText(index)={['MMSEaproxPCG']}; end
        if flagFMMSE, semilogy(turboAxis,[BER_soft_FMMSE(snr),BER_soft_FMMSE_turbo(:,snr).'],'o-m'),hold on,...
                index=index+1; legendText(index)={['FMMSE']}; end
        if flagBEP, semilogy(turboAxis,[BER_soft_BEP(snr),BER_soft_BEP_turbo(:,snr).'],'s-r'),hold on,...
                index=index+1; legendText(index)={['BEP']}; end
        if flagPBEP, semilogy(turboAxis,[BER_soft_PBEP(snr),BER_soft_PBEP_turbo(:,snr).'],'v-r'),hold on,...
                index=index+1; legendText(index)={['P-BEP']}; end
        if flagDBEP, semilogy(turboAxis,[BER_soft_DBEP(snr),BER_soft_DBEP_turbo(:,snr).'],'d-r'),hold on,...
                index=index+1; legendText(index)={['D-BEP']}; end
        if flagDBEPAproxGS, semilogy(turboAxis,[BER_soft_DBEPAproxGS(snr),BER_soft_DBEPAproxGS_turbo(:,snr).'],'*-c'),hold on,...
                index=index+1; legendText(index)={['D-BEPaproxGS']}; end
        if flagDBEPAproxNeuman, semilogy(turboAxis,[BER_soft_DBEPAproxNeuman(snr),BER_soft_DBEPAproxNeuman_turbo(:,snr).'],'*-k'),hold on,...
                index=index+1; legendText(index)={['D-BEPaproxNeuman']}; end
        if flagDBEPAproxPCG, semilogy(turboAxis,[BER_soft_DBEPAproxPCG(snr),BER_soft_DBEPAproxPCG_turbo(:,snr).'],'*-g'),hold on,...
                index=index+1; legendText(index)={['D-BEPaproxPCG']}; end
        if flagPFEP, semilogy(turboAxis,[BER_soft_PFEP(snr),BER_soft_PFEP_turbo(:,snr).'],'d-c'),hold on,...
                index=index+1; legendText(index)={['P-FEP']}; end
        if flagDFEP, semilogy(turboAxis,[BER_soft_DFEP(snr),BER_soft_DFEP_turbo(:,snr).'],'*-c'),hold on,...
                index=index+1; legendText(index)={['D-FEP']}; end
        if flagSEP, semilogy(turboAxis,[BER_soft_SEP(snr),BER_soft_SEP_turbo(:,snr).'],'o-g'),hold on,...
                index=index+1; legendText(index)={['SEP']}; end
        if flagPKSEP, semilogy(turboAxis,[BER_soft_PKSEP(snr),BER_soft_PKSEP_turbo(:,snr).'],'*-g'),hold on,...
                index=index+1; legendText(index)={['P-KSEP']}; end
        if flagDKSEP, semilogy(turboAxis,[BER_soft_DKSEP(snr),BER_soft_DKSEP_turbo(:,snr).'],'*-r'),hold on,...
                index=index+1; legendText(index)={['D-KSEP']}; end
        if flagBEPSun, semilogy(turboAxis,[BER_soft_BPEP(snr),BER_soft_BPEP_turbo(:,snr).'],'+-c'),hold on,...
                index=index+1; legendText(index)={['BPEP']}; end
        if flagEPICLMMSE, semilogy(turboAxis,[BER_soft_EPICLMMSE(snr),BER_soft_EPICLMMSE_turbo(:,snr).'],'+-m'),hold on,...
                index=index+1; legendText(index)={['EP IC-LMMSE']}; end
        if flagOptimal, semilogy(turboAxis,[BER_soft_Optimal(snr),BER_soft_Optimal_turbo(:,snr).'],'d-k'),hold on,...
                index=index+1; legendText(index)={['Optimal']}; end
        
        hold off
        legend(legendText)
        xlabel('turbo iteration (t)'),...
            ylabel('BER')
        
        if dataEP.scenario==1
            title([num2str(dataEP.numberChannels),' channels with ',num2str(dataEP.numberTaps),...
                ' taps, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
                ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
                ' frames, ', constelation_name, ' symbols and SNRdB=', num2str(SNRdBSoft(snr)), 'dB']);
        else
            title([num2str(dataEP.numberChannels),' channels with Ntx=',num2str(dataEP.numberAntennas(1)),...
                ' and Nrx=',num2str(dataEP.numberAntennas(2)),...
                ' antennas, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
                ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
                ' frames and ', constelation_name, ' symbols and SNRdB=', num2str(SNRdBSoft(snr)), 'dB']);
        end
        
        axis([turboAxis(1) turboAxis(end) 1e-6 1])
        
    end
    
    
end

if dataEP.flagIterDec
    turboAxis=0:dataEP.numberTurbo;
    for snr=1:length(SNRdBSoft)
        figure,
        index=0;
        clearvars legendText
        
        if flagMMSE, plot(turboAxis,niterdec_MMSE(:,snr),'v-b'),hold on,...
                index=index+1; legendText(index)={['MMSE']};  end
        if flagMMSEAproxGS, semilogy(turboAxis,niterdec_MMSEAproxGS(:,snr),'o-g'),hold on,...
                index=index+1; legendText(index)={['MMSEaproxGS']}; end
        if flagMMSEAproxNeuman, semilogy(turboAxis,niterdec_MMSEAproxNeuman(:,snr),'o-m'),hold on,...
                index=index+1; legendText(index)={['MMSEaproxNeuman']}; end
        if flagMMSEAproxPCG, semilogy(turboAxis,niterdec_MMSEAproxPCG(:,snr),'o-r'),hold on,...
                index=index+1; legendText(index)={['MMSEaproxPCG']}; end
        if flagFMMSE, semilogy(turboAxis,niterdec_FMMSE(:,snr),'o-m'),hold on,...
                index=index+1; legendText(index)={['FMMSE']}; end
        if flagBEP, semilogy(turboAxis,niterdec_BEP(:,snr),'s-r'),hold on,...
                index=index+1; legendText(index)={['BEP']}; end
        if flagPBEP, semilogy(turboAxis,niterdec_PBEP(:,snr),'v-r'),hold on,...
                index=index+1; legendText(index)={['P-BEP']}; end
        if flagDBEP, semilogy(turboAxis,niterdec_DBEP(:,snr),'d-r'),hold on,...
                index=index+1; legendText(index)={['D-BEP']}; end
        if flagDBEPAproxGS, semilogy(turboAxis,niterdec_DBEPAproxGS(:,snr),'*-c'),hold on,...
                index=index+1; legendText(index)={['D-BEPaproxGS']}; end
        if flagDBEPAproxNeuman, semilogy(turboAxis,niterdec_DBEPAproxNeuman(:,snr),'*-k'),hold on,...
                index=index+1; legendText(index)={['D-BEPaproxNeuman']}; end
        if flagDBEPAproxPCG, semilogy(turboAxis,niterdec_DBEPAproxPCG(:,snr),'*-g'),hold on,...
                index=index+1; legendText(index)={['D-BEPaproxPCG']}; end
        if flagPFEP, semilogy(turboAxis,niterdec_PFEP(:,snr),'d-c'),hold on,...
                index=index+1; legendText(index)={['P-FEP']}; end
        if flagDFEP, semilogy(turboAxis,niterdec_DFEP(:,snr),'*-c'),hold on,...
                index=index+1; legendText(index)={['D-FEP']}; end
        if flagSEP, semilogy(turboAxis,niterdec_SEP(:,snr),'o-g'),hold on,...
                index=index+1; legendText(index)={['SEP']}; end
        if flagPKSEP, semilogy(turboAxis,niterdec_PKSEP(:,snr),'*-g'),hold on,...
                index=index+1; legendText(index)={['P-KSEP']}; end
        if flagDKSEP, semilogy(turboAxis,niterdec_DKSEP(:,snr),'*-r'),hold on,...
                index=index+1; legendText(index)={['D-KSEP']}; end
        if flagBEPSun, semilogy(turboAxis,niterdec_BPEP(:,snr),'+-c'),hold on,...
                index=index+1; legendText(index)={['BPEP']}; end
        if flagEPICLMMSE, semilogy(turboAxis,niterdec_EPICLMMSE(:,snr),'+-m'),hold on,...
                index=index+1; legendText(index)={['EP IC-LMMSE']}; end
        if flagOptimal, semilogy(turboAxis,niterdec_Optimal(:,snr),'d-k'),hold on,...
                index=index+1; legendText(index)={['Optimal']}; end
        
        hold off
        legend(legendText)
        xlabel('turbo iteration (t)'),...
            ylabel('number of iterations of the LDPC')
        
        if dataEP.scenario==1
            title([num2str(dataEP.numberChannels),' channels with ',num2str(dataEP.numberTaps),...
                ' taps, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
                ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
                ' frames, ', constelation_name, ' symbols and SNRdB=', num2str(SNRdBSoft(snr)), 'dB']);
        else
            title([num2str(dataEP.numberChannels),' channels with Ntx=',num2str(dataEP.numberAntennas(1)),...
                ' and Nrx=',num2str(dataEP.numberAntennas(2)),...
                ' antennas, entr cod=',num2str(2^dataEP.bits),', sal cod=',num2str(2^dataEP.bits/rate),...
                ', ', num2str(dataEP.numberFrames*dataEP.numberSimulations),...
                ' frames and ', constelation_name, ' symbols and SNRdB=', num2str(SNRdBSoft(snr)), 'dB']);
        end
    end
end

if dataEP.numberTurbo>0
save(saveName,...
    'SNRdBSoft','dataEP',...
    'BER_soft_BEP','BER_soft_BEP_turbo',...
    'BER_soft_PBEP','BER_soft_PBEP_turbo',...
    'BER_soft_DBEP','BER_soft_DBEP_turbo',...
    'BER_soft_MMSE','BER_soft_MMSE_turbo',...
    'BER_soft_FMMSE','BER_soft_FMMSE_turbo',...
    'BER_soft_PFEP','BER_soft_PFEP_turbo',...
    'BER_soft_DFEP','BER_soft_DFEP_turbo',...
    'BER_soft_BPEP','BER_soft_BPEP_turbo',...
    'BER_soft_SEP','BER_soft_SEP_turbo',...
    'BER_soft_MMSEAproxGS','BER_soft_MMSEAproxGS_turbo',...
    'BER_soft_MMSEAproxNeuman','BER_soft_MMSEAproxNeuman_turbo',...
    'BER_soft_DBEPAproxNeuman','BER_soft_DBEPAproxNeuman_turbo',...
    'BER_soft_MMSEAproxPCG','BER_soft_MMSEAproxPCG_turbo',...
    'BER_soft_DBEPAproxPCG','BER_soft_DBEPAproxPCG_turbo',...
    'BER_soft_DBEPAproxGS','BER_soft_DBEPAproxGS_turbo',...    
    'BER_soft_DKSEP','BER_soft_DKSEP_turbo',...
    'BER_soft_EPICLMMSE','BER_soft_EPICLMMSE_turbo',...
    'BER_soft_PKSEP','BER_soft_PKSEP_turbo',...
    'BER_hard_PBEP','BER_hard_PBEP_turbo',...
    'BER_hard_DBEP','BER_hard_DBEP_turbo',...
    'BER_hard_MMSE','BER_hard_MMSE_turbo',...
    'BER_hard_FMMSE','BER_hard_FMMSE_turbo',...
    'BER_hard_SEP','BER_hard_SEP_turbo',...
    'BER_hard_MMSEAproxGS','BER_hard_MMSEAproxGS_turbo',...
    'BER_hard_MMSEAproxNeuman','BER_hard_MMSEAproxNeuman_turbo',...
    'BER_hard_MMSEAproxPCG','BER_hard_MMSEAproxPCG_turbo',...
    'BER_hard_EPICLMMSE','BER_hard_EPICLMMSE_turbo')
%'BER_hard_PFEP','BER_hard_PFEP_turbo',...
%'BER_hard_DFEP','BER_hard_DFEP_turbo',...
%'BER_hard_BPEP','BER_hard_BPEP_turbo',...
%'BER_hard_DKSEP','BER_hard_DKSEP_turbo',...
%'BER_hard_DBEPAproxPCG','BER_hard_DBEPAproxPCG_turbo',...
%'BER_hard_DBEPAproxNeuman','BER_hard_DBEPAproxNeuman_turbo',...
% 'BER_hard_PKSEP','BER_hard_PKSEP_turbo'

else
    
save(saveName,...
    'SNRdBSoft','dataEP',...
    'BER_soft_BEP',...    
    'BER_soft_PBEP',...
    'BER_soft_DBEP',...
    'BER_soft_MMSE',...
    'BER_soft_FMMSE',...
    'BER_soft_PFEP',...
    'BER_soft_DFEP',...
    'BER_soft_BPEP',...
    'BER_soft_SEP',...
    'BER_soft_MMSEAproxGS',...
    'BER_soft_MMSEAproxNeuman',...
    'BER_soft_DBEPAproxNeuman',...
    'BER_soft_MMSEAproxPCG',...
    'BER_soft_DBEPAproxPCG',...
    'BER_soft_DBEPAproxGS',...    
    'BER_soft_DKSEP',...
    'BER_soft_EPICLMMSE',...
    'BER_soft_PKSEP',...
    'BER_hard_BEP',...
    'BER_hard_PBEP',...
    'BER_hard_DBEP',...
    'BER_hard_MMSE',...
    'BER_hard_FMMSE',...
    'BER_hard_SEP',...
    'BER_hard_MMSEAproxGS',...
    'BER_hard_MMSEAproxNeuman',...
    'BER_hard_MMSEAproxPCG',...
    'BER_hard_DBEPAproxPCG',...
    'BER_hard_DBEPAproxGS',...    
    'BER_hard_DKSEP',...    
    'BER_hard_EPICLMMSE',...
    'BER_hard_PKSEP')
%'BER_hard_PFEP','BER_hard_PFEP_turbo',...
%'BER_hard_DFEP','BER_hard_DFEP_turbo',...
%'BER_hard_BPEP','BER_hard_BPEP_turbo',...
%'BER_hard_DKSEP','BER_hard_DKSEP_turbo',...
%'BER_hard_DBEPAproxPCG','BER_hard_DBEPAproxPCG_turbo',...
%'BER_hard_DBEPAproxNeuman','BER_hard_DBEPAproxNeuman_turbo',...
% 'BER_hard_PKSEP','BER_hard_PKSEP_turbo'
end
