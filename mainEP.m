function mainEP(flagRun,dataEP,scenario)

if nargin==0
    flagRun=1;
end
if dataEP.flagPSK, modName='PSK'; else modName='QAM'; end

if dataEP.scenario==1
    dataEP.folderSave=['ResultsEQ/Results'];
else %MIMO
    dataEP.folderSave=['ResultsMIMO/Results'];
end

if dataEP.scenario==2, Ants=['MIMO',num2str(dataEP.numberAntennas(1)),'x',num2str(dataEP.numberAntennas(2))],  end


dataEP.folderSave=[dataEP.folderSave,num2str(dataEP.M),modName,'Lim',num2str(dataEP.LLRlim),dataEP.channelName];

if dataEP.flagBlockLengths==1
    dataEP.folderSave=[dataEP.folderSave,'N',num2str(dataEP.channelBlockLength)];
end

if dataEP.flagEPiterations
    dataEP.folderSave=[dataEP.folderSave,'S',num2str(dataEP.EPiterations)];
    dataEP.BEP_S=dataEP.EPiterations;% BEP [Santos16], Default 10
    dataEP.PBEP_S=dataEP.EPiterations; % P-BEP [Santos18], Default 3
    dataEP.DBEP_S=dataEP.EPiterations; % D-BEP [Santos19], Default 1
    dataEP.KSEP_S=dataEP.EPiterations;
%     dataEP.DKSEP_S=dataEP.EPiterations;
% else
%     dataEP.BEP_S=10;  % BEP [Santos16], Default 10
%     dataEP.PBEP_S=3;  % P-BEP [Santos18], Default 3
%     dataEP.DBEP_S=1;  % D-BEP [Santos19], Default 1
%     dataEP.KSEP_S=3;
%     dataEP.DKSEP_S=1; % D-KSEP [Santos19], Default 1
end

% 
%% Executing the simulation
if flagRun
 parfor simulationNumber=1:dataEP.numberSimulations
     %dataEP.simulationNumber=simulationNumber;
     runSimulationNumber(dataEP,scenario,simulationNumber);
 end
end

%% Plotting the figures
nameSave=[dataEP.folderSave,'/BER'];


if dataEP.flagEPiterations
nameSave=[nameSave,'vsTS',num2str(dataEP.EPiterations)];
else
nameSave=[nameSave,'vsSNRdB'];    
end

if dataEP.flagBlockLengths
nameSave=[nameSave,'N',num2str(dataEP.channelBlockLength)];
end

if dataEP.scenario==1
    Ants=[];
end
if isfield(dataEP,'long'),
    aux=num2str(dataEP.rate,2);Rate=aux(3:end);
    nameSave=[nameSave,num2str(dataEP.M),modName,Ants,'Lim',num2str(dataEP.LLRlim),'Long',Rate,dataEP.channelName,'Fr',num2str(dataEP.numberFrames*dataEP.numberSimulations)];
else  
    nameSave=[nameSave,num2str(dataEP.M),modName,Ants,'Lim',num2str(dataEP.LLRlim),dataEP.channelName,'Fr',num2str(dataEP.numberFrames*dataEP.numberSimulations)];
end


plotEPeqSNRdBSave(dataEP,nameSave,scenario)



