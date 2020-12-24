function runSimulationNumber(dataEP,scenario,simulationNumber)
%
% Function: runSimulationNumber(dataEP,scenario,simulationNumber)
%
%
% Author: Irene Santos, Juan Jose Murillo-Fuentes
%
% Contact: murillo@us.es, irenesantos@us.es
%
% Created 20/12/2018
%
% Description: This function runs different equalizers (MMSE, EP,
% Optimal) and outputs the results into a text file.
%
% Inputs:
%   dataEP: is a struct with many parameters defining the simulation,
%         please, see the configuration files. Some fields are:
%     dataEP.numberFrames is the number of words transmitted per simulation
%     dataEP.bits is the size of information dataEP.bits
%     dataEP.numberTaps is the number of taps of the channel in equalization
%     dataEP.numberAntennas is the number of transmitted and receive antennas in MIMO
%     dataEP.numberChannels is the number of simulated channels
%     dataEP.M is the constellation size
%     dataEP.flagPSK is 1 if the complex constellation is a PSK and 0 if it is a QAM
%     dataEP.flagIterDec is 1 if the figures of number of iterations of the channel
%     decoder is plotted
%   scenario: determines if MIMO or equalization is run, the range of
%     dB to simulate, please, see the configuration files
%   simulationNumber:  allows for several simulations run in paralel, 
%     indicating the number of the simulation. 

dataEP.simulationNumber=simulationNumber;

%% NAME OF OUTPUT FILE to save detailed results
%To write and later read results
DELIMITER=' ';
HEADERLINES=4;

%JJMF meti un folder
if ~isfield(dataEP,'folderSave'), dataEP.folderSave='Results'; end
[status, msg, msgID] = mkdir(dataEP.folderSave);

%JJMF No entiendo lo que sale por pantalla. ?Por que cuando termina con el
%canal 10/10 vuelve a empezar? ?Cuantas veces lo repite?

%JJMF Puse un smart indent

%JJMF Me pregunto si hay posibilidad de seleccionar qu? m?todos se simulan

if dataEP.scenario==1 % equalization
    nameO=([dataEP.folderSave,'/','resultsEPeq', '_b' num2str(dataEP.bits) '_M' num2str(dataEP.M) '_L' num2str(dataEP.numberTaps) '_n' num2str(dataEP.channelBlockLength) '_numCh' num2str(dataEP.numberChannels) '_' num2str(dataEP.simulationNumber) '.txt']);
    nameMSE=([dataEP.folderSave,'/','MSEresultsEPeq', '_b' num2str(dataEP.bits) '_M' num2str(dataEP.M) '_L' num2str(dataEP.numberTaps) '_n' num2str(dataEP.channelBlockLength) '_numCh' num2str(dataEP.numberChannels) '_' num2str(dataEP.simulationNumber) '.txt']);
else
    nameO=([dataEP.folderSave,'/','resultsEPmm', '_b' num2str(dataEP.bits) '_M' num2str(dataEP.M) '_NtxNr' num2str(dataEP.numberAntennas(1)) 'x' num2str(dataEP.numberAntennas(2)) '_n' num2str(dataEP.channelBlockLength) '_numCh' num2str(dataEP.numberChannels) '_' num2str(dataEP.simulationNumber) '.txt']);
    nameMSE=([dataEP.folderSave,'/','MSEresultsEPmm', '_b' num2str(dataEP.bits) '_M' num2str(dataEP.M) '_NtxNr' num2str(dataEP.numberAntennas(1)) 'x' num2str(dataEP.numberAntennas(2)) '_n' num2str(dataEP.channelBlockLength) '_numCh' num2str(dataEP.numberChannels) '_' num2str(dataEP.simulationNumber) '.txt']);
end

%% Including the subfolder with the code of the algorihtms to the path
algstuffroot=[pwd '/algorithmsCode/'];
addpath(genpath(algstuffroot))




%% INITIALIZATION WITH INPUT PARAMETERS

verbose=1; % To display some text when simulating

if dataEP.scenario==1 % equalization
    nameLog=([dataEP.folderSave,'/','LogEPeq', '_b' num2str(dataEP.bits) '_L' num2str(dataEP.numberTaps) '_numCh' num2str(dataEP.numberChannels) '_' num2str(dataEP.simulationNumber) '.txt']);
else
    nameLog=([dataEP.folderSave,'/','LogEPmm', '_b' num2str(dataEP.bits) '_M' num2str(dataEP.M) '_NtxNr' num2str(dataEP.numberAntennas(1)) 'x' num2str(dataEP.numberAntennas(2)) '_numCh' num2str(dataEP.numberChannels) '_' num2str(dataEP.simulationNumber) '.txt']);
end

%
if verbose
    fidLog=fopen(nameLog,'w');
end

%dataEP.numberFrames
% Parameter dataEP.numberFrames
% We send dataEP.numberFrames frames
if verbose, disp('Setting numberFrames');
    count=fwrite(fidLog,'Setting numberFrames','char');
    count=fwrite(fidLog,10,'uint8');
end

%bits
%parameter dataEP.bits, to choose the parity matrix, it is the number of bits of
%the transmitted word. Rate 1/2 is used, with LDPC codes
if isfield(dataEP,'long')
    H = dvbs2ldpc(dataEP.rate);
    %H=full(H);
    rate=dataEP.rate;
else
    H=loadParityCheckMatrix(2^dataEP.bits);
    rate=1/2;
if verbose, disp('Parity check matrix file loaded');
    count=fwrite(fidLog,'Parity check matrix file loaded','char');
    count=fwrite(fidLog,10,'uint8');
    
end
end
encoderOutputLenght=size(H,2);

encoderInputLenght=encoderOutputLenght*rate;

%dataEP.M
%Parameter dataEP.M is the size of the constellation. Must be a multiple of 2
Mnew=2^ceil(log2(dataEP.M));
if Mnew~=dataEP.M
    if verbose, disp('M rounded to upper multiple of 2');
        count=fwrite(fidLog,'M rounded to upper multiple of 2','char');
        count=fwrite(fidLog,10,'uint8');
    end
    dataEP.M=Mnew;
end

%With dataEP.M we initializate the following
%Modulation and coding -> size of frames(=words)
modulatorOutputLength=encoderOutputLenght/log2(dataEP.M);
modulatorOutputLengthNew=log2(dataEP.M)*ceil(modulatorOutputLength/log2(dataEP.M)); % we round modulatorOutputLength to the nearest integer
encoderOutputLenght_new=modulatorOutputLengthNew*log2(dataEP.M);

% Parameter channelBlockLength is the size of channel in equalization (we divide the whole word in words of size channelBlockLength)
% We do not use it in MIMO
if dataEP.channelBlockLength>modulatorOutputLengthNew
    dataEP.channelBlockLength=modulatorOutputLengthNew;
end
if mod(modulatorOutputLengthNew,dataEP.channelBlockLength)~=0
    dataEP.channelBlockLength=modulatorOutputLengthNew/round(modulatorOutputLengthNew/dataEP.channelBlockLength);
    if verbose
        disp('dimCh rounded to the nearest multiple of modulatorOutputLengthNew'),
        count=fwrite(fidLog,'dimCh rounded to be multiple of modulatorOutputLengthNew','char');
        count=fwrite(fidLog,10,'uint8');
    end
end

%dataEP.numberTaps or numberAntennas
%Parameter numberTaps is the length of the channel (in
%equalization). dataEP.numberTaps=1 if no
%memory or SISO system.
% numberAntennas is the number of antennas (in MIMO).
if dataEP.numberTaps<0 || dataEP.numberAntennas(1)<0 || dataEP.numberAntennas(2)<0
    if verbose, disp('L set to two'),
        count=fwrite(fidLog,'L set to two','char');
        count=fwrite(fidLog,10,'uint8');
    end
    dataEP.numberTaps=2;
    dataEP.numberAntennas=[2 2];
end

% SNRdBIni, SNRdBStep, SNRdBEnd
% To define the snr range to simulate
SNRdB=dataEP.SNRdBIni:dataEP.SNRdBStep:dataEP.SNRdBEnd;

%Scenario. 
if dataEP.scenario==1 % equalization case
    flagMIMO=0; % For equalization this flag is zero
else % MIMO case
    flagMIMO=1;
end

% Parameter numberAntennas(1) is the size of each sent word in MIMO
% (we divide the whole word in words of size numberAntennas(1))
if flagMIMO
    modulatorOutputLengthNew=lcm(log2(dataEP.M),dataEP.numberAntennas(1))*ceil(modulatorOutputLength/lcm(log2(dataEP.M),dataEP.numberAntennas(1)));
    encoderOutputLenght_new=modulatorOutputLengthNew*log2(dataEP.M);
end

% dataEP.channel_h
% If channel_h=[] we average over numberChannels different channels. If a
% value is given to channel_h, the number of rows is the number of different
% channels and the number of columns is the number of taps.
if ~isempty(dataEP.channel_h) && ~flagMIMO
    [dataEP.numberChannels,numberTaps]=size(dataEP.channel_h);
    if verbose
        disp('No random channels'),
        count=fwrite(fidLog,'No random channels','char');
        count=fwrite(fidLog,10,'uint8');
    end
end

% numberChannel
% numerChannel is the number of the random channels simulated. It is used
% for the seed of the random number generator
idxChannels=1:dataEP.numberChannels;

%dataEP.numberTurbo
% If numberTurbo<0, we rounded it to zero.
if dataEP.numberTurbo<0
    dataEP.numberTurbo=0;
end


%% OTHER INITIALIZATION 
if ~isfield(dataEP,'LLRlim'),  limLDPC=3;else
    limLDPC=dataEP.LLRlim; end
    % We limit the LLR given to the LDPC

% LDPC encoder
hEnc = comm.LDPCEncoder('ParityCheckMatrix',sparse(H));
% Modulator
if dataEP.complexFlag==0
    hMod = comm.PAMModulator(dataEP.M, 'BitInput',true,'SymbolMapping','Gray');
    hDemod = comm.PAMDemodulator(dataEP.M, 'BitOutput',true,'SymbolMapping','Gray');
else
    if dataEP.flagPSK
        hMod = comm.PSKModulator('ModulationOrder',dataEP.M,'BitInput',true,'SymbolMapping','Gray');
        hDemod = comm.PSKDemodulator('ModulationOrder',dataEP.M,'BitOutput',true,'SymbolMapping','Gray');
    else
        hMod = comm.RectangularQAMModulator('ModulationOrder',dataEP.M,'BitInput',true,'SymbolMapping','Gray');
        hDemod = comm.RectangularQAMDemodulator('ModulationOrder',dataEP.M,'BitOutput',true,'SymbolMapping','Gray');
    end
end

% Interleaving/deinterleaving
seed_inter=4831; % random stream that determines the specific permutation

% LDPC decoder
hDec = comm.LDPCDecoder('ParityCheckMatrix',hEnc.ParityCheckMatrix,'DecisionMethod','Soft decision','OutputValue','Whole codeword','MaximumIterationCount',dataEP.maxIterLDPC,'IterationTerminationCondition','Parity check satisfied','NumIterationsOutputPort',true);


%% Generation of the constellation
% We use the comm toolbox, where the points of the constellation are
% separated by two
A=transpose(hMod.constellation);
energy=sum(abs(A).^2)/dataEP.M; %Average Symbol Energy
%dd=1;Acons=dd^2*(2*dataEP.M-2)/3

symbolmap=reshape(step(hDemod,A.'),log2(dataEP.M),dataEP.M); % Column i is the gray representation of symbol A(i)


%% Computing the variance of the noise from the SNRdB (abcisa axis)

if flagMIMO % MIMO, SNRdB is treated as EsNo
    esno=10.^(SNRdB/10);%*log2(dataEP.M); %ebno in the MIMO case is snr
    
    % Sigma is different for hard and soft. We simulate with sigma without
    % rate, i.e., sigma for hard estimations. After simulation, we should
    % move the soft curve 3dB to the right.
    sigma=sqrt(dataEP.numberAntennas(1)*energy./(esno*2));
else % Equalization
    ebno=10.^(SNRdB/10);
    esno=log2(dataEP.M)*ebno;
    EsNo=10*log10(esno);
    sigma=10.^(-EsNo/20)*sqrt(energy/2);
end

% If error in the channel matrix we use a modified noise variance
sigma_novarCH=sigma;
sigma=sqrt(sigma.^2+dataEP.numberAntennas(1)*dataEP.varnoiseCH*energy);

%%%

%% Creating new file
% and writting data of the scenario
fidauxMSE=fopen(nameMSE,'w');
count=fwrite(fidauxMSE,'numberFrames=','char');
count=fwrite(fidauxMSE,num2str(dataEP.numberFrames),'char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,'b=','char');
count=fwrite(fidauxMSE,num2str(dataEP.bits),'char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,'numerChannels=','char');
count=fwrite(fidauxMSE,num2str(dataEP.numberChannels),'char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,'M=','char');
count=fwrite(fidauxMSE,num2str(dataEP.M),'char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,'complexFlag=','char');
count=fwrite(fidauxMSE,num2str(dataEP.complexFlag),'char');
count=fwrite(fidauxMSE,' ','char');
if dataEP.scenario==1
    count=fwrite(fidauxMSE,'L=','char');
    count=fwrite(fidauxMSE,num2str(dataEP.numberTaps),'char');
else
    count=fwrite(fidauxMSE,'Nt=','char');
    count=fwrite(fidauxMSE,num2str(dataEP.numberAntennas(1)),'char');
    count=fwrite(fidauxMSE,' ','char');
    count=fwrite(fidauxMSE,'Nr=','char');
    count=fwrite(fidauxMSE,num2str(dataEP.numberAntennas(2)),'char');
end
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,'scenario=','char');
count=fwrite(fidauxMSE,num2str(dataEP.scenario),'char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,'SNRdB=','char');
count=fwrite(fidauxMSE,num2str(dataEP.SNRdBIni),'char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,num2str(dataEP.SNRdBStep),'char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,num2str(dataEP.SNRdBEnd),'char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,'numberTurbo=','char');
count=fwrite(fidauxMSE,num2str(dataEP.numberTurbo),'char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,10,'uint8');

count=fwrite(fidauxMSE,'SNRdB=SNRdBIni:step:SNRdBEnd','char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,'numberTurbo','char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,'numberChannels','char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,'numberTaps','char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,10,'uint8');

count=fwrite(fidauxMSE,num2str(dataEP.SNRdBIni),'char');
count=fwrite(fidauxMSE,' ','uint8');
count=fwrite(fidauxMSE,num2str(dataEP.SNRdBStep),'char');
count=fwrite(fidauxMSE,' ','uint8');
count=fwrite(fidauxMSE,num2str(dataEP.SNRdBEnd),'char');
count=fwrite(fidauxMSE,' ','uint8');
count=fwrite(fidauxMSE,num2str(dataEP.numberTurbo),'char');
count=fwrite(fidauxMSE,' ','uint8');
count=fwrite(fidauxMSE,num2str(dataEP.numberChannels),'char');
count=fwrite(fidauxMSE,' ','uint8');
if dataEP.scenario==1
    count=fwrite(fidauxMSE,num2str(dataEP.numberTaps),'char');
else
    count=fwrite(fidauxMSE,num2str(dataEP.numberAntennas(1)),'char');
    count=fwrite(fidauxMSE,' ','uint8');
    count=fwrite(fidauxMSE,num2str(dataEP.numberAntennas(2)),'char');
end
count=fwrite(fidauxMSE,' ','uint8');
count=fwrite(fidauxMSE,10,'uint8');

count=fwrite(fidauxMSE,'n_h','char');
count=fwrite(fidauxMSE,' ','char');
count=fwrite(fidauxMSE,'hDBEPC','char');
count=fwrite(fidauxMSE,' ','uint8');
count=fwrite(fidauxMSE,'hDBEPD','char');
count=fwrite(fidauxMSE,' ','uint8');
count=fwrite(fidauxMSE,'hDBEPL','char');
count=fwrite(fidauxMSE,' ','uint8');
count=fwrite(fidauxMSE,10,'uint8');
fclose(fidauxMSE);
%%%







%% Creating new file
% and writting data of the scenario
fidaux=fopen(nameO,'w');
count=fwrite(fidaux,'numberFrames=','char');
count=fwrite(fidaux,num2str(dataEP.numberFrames),'char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,'b=','char');
count=fwrite(fidaux,num2str(dataEP.bits),'char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,'numerChannels=','char');
count=fwrite(fidaux,num2str(dataEP.numberChannels),'char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,'M=','char');
count=fwrite(fidaux,num2str(dataEP.M),'char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,'complexFlag=','char');
count=fwrite(fidaux,num2str(dataEP.complexFlag),'char');
count=fwrite(fidaux,' ','char');
if dataEP.scenario==1
    count=fwrite(fidaux,'L=','char');
    count=fwrite(fidaux,num2str(dataEP.numberTaps),'char');
else
    count=fwrite(fidaux,'Nt=','char');
    count=fwrite(fidaux,num2str(dataEP.numberAntennas(1)),'char');
    count=fwrite(fidaux,' ','char');
    count=fwrite(fidaux,'Nr=','char');
    count=fwrite(fidaux,num2str(dataEP.numberAntennas(2)),'char');
end
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,'scenario=','char');
count=fwrite(fidaux,num2str(dataEP.scenario),'char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,'SNRdB=','char');
count=fwrite(fidaux,num2str(dataEP.SNRdBIni),'char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,num2str(dataEP.SNRdBStep),'char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,num2str(dataEP.SNRdBEnd),'char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,'numberTurbo=','char');
count=fwrite(fidaux,num2str(dataEP.numberTurbo),'char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,10,'uint8');

count=fwrite(fidaux,'SNRdB=SNRdBIni:step:SNRdBEnd','char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,'numberTurbo','char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,'numberChannels','char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,'numberTaps','char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,10,'uint8');

count=fwrite(fidaux,num2str(dataEP.SNRdBIni),'char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,num2str(dataEP.SNRdBStep),'char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,num2str(dataEP.SNRdBEnd),'char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,num2str(dataEP.numberTurbo),'char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,num2str(dataEP.numberChannels),'char');
count=fwrite(fidaux,' ','uint8');
if dataEP.scenario==1
    count=fwrite(fidaux,num2str(dataEP.numberTaps),'char');
else
    count=fwrite(fidaux,num2str(dataEP.numberAntennas(1)),'char');
    count=fwrite(fidaux,' ','uint8');
    count=fwrite(fidaux,num2str(dataEP.numberAntennas(2)),'char');
end
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,10,'uint8');

count=fwrite(fidaux,'n_h','char');
count=fwrite(fidaux,' ','char');
count=fwrite(fidaux,'hMMSE','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hMMSEturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hMMSEAproxGS','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hMMSEAproxGSturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hMMSEAproxNeuman','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hMMSEAproxNeumanturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hMMSEAproxPCG','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hMMSEAproxPCGturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hFMMSE','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hFMMSEturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hBEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hBEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hPBEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hPBEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hDBEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hDBEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hDBEPAprox','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hDBEPAproxturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hPFEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hPFEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hDFEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hDFEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hSEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hSEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hPKSEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hPKSEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hEPICLMMSE','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hEPICLMMSEturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hOptimal','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'hOptimalturbo','char');
count=fwrite(fidaux,' ','uint8');

count=fwrite(fidaux,'sMMSE','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sMMSEturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sMMSEAproxGS','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sMMSEAproxGSturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sMMSEAproxNeuman','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sMMSEAproxNeumanturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sMMSEAproxPCG','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sMMSEAproxPCGturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sFMMSE','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sFMMSEturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sBEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sBEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sPBEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sPBEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sDBEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sDBEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sDBEPAprox','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sDBEPAproxturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sPFEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sPFEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sDFEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sDFEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sSEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sSEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sPKSEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sPKSEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sDKSEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sKDSEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sEPICLMMSE','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sEPICLMMSEturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sBPEP','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sBPEPturbo','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sOptimal','char');
count=fwrite(fidaux,' ','uint8');
count=fwrite(fidaux,'sOptimalturbo','char');
count=fwrite(fidaux,10,'uint8');
fclose(fidaux);


%% RUNNING MONTECARLO
tStart=tic;

for n_h=idxChannels
    
    %JJMF         
    if verbose, disp(['SimulationNumber =',num2str(dataEP.simulationNumber)])%),...
          %  ', n_h=',num2str(n_h),'/',num2str(dataEP.numberChannels)])
          end
            
    %%Computing a sample of the channel
    if flagMIMO % MIMO
       
        % Seed
       % randn('seed',n_h);
       % rand('seed',n_h);
        rng(n_h,'twister')
        
        if dataEP.complexFlag==0 % real
            % Gaussian channels
            Hmat_novarCH=randn(dataEP.numberAntennas(2),dataEP.numberAntennas(1));
            Hmat=Hmat_novarCH+randn(size(Hmat_novarCH))*sqrt(dataEP.varnoiseCH);
        else
            Hmat_novarCH=randn(dataEP.numberAntennas(2),dataEP.numberAntennas(1))+1i*randn(dataEP.numberAntennas(2),dataEP.numberAntennas(1));
            Hmat_novarCH=Hmat_novarCH/sqrt(2);
            Hmat=Hmat_novarCH+randn(size(Hmat_novarCH))*sqrt(dataEP.varnoiseCH)+1i*randn(size(Hmat_novarCH))*sqrt(dataEP.varnoiseCH);
        end
%         delta=0.2;Hmat=[ones(1,4);zeros(3,4)];Hmat(:,2)=Hmat(:,2)+delta;Hmat(:,3)=Hmat(:,3)+2*delta;Hmat(:,4)=Hmat(:,4)+4*delta;
%         Hmat_novarCH=Hmat;
        %Hmat_novarCH=[1;1];
%         Hmat_novarCH=ones(32,32)*0.5+.5*eye(32);
%         Hmat=Hmat_novarCH;
%         Hmat
    else %EQUALIZATION
        
        % Channel
        if isempty(dataEP.channel_h)
            % Seed
            %randn('state',100*n_h*7);
            %rand('state',100*n_h*7);
            rng(100*n_h*7,'twister')
            
            if dataEP.complexFlag==0 % real
                % Gaussian channels
                h_novarCH=randn(1,dataEP.numberTaps)/sqrt(dataEP.numberTaps);
                if scenario.flagNorm  %Normalized if set in the scenario
                    h_novarCH=h_novarCH./norm(h_novarCH);
                end
                h=h_novarCH+randn(size(h_novarCH))*sqrt(dataEP.varnoiseCH);
            else
                % Gaussian channels
                h_novarCH=randn(1,dataEP.numberTaps)+1i*randn(1,dataEP.numberTaps);
                h_novarCH=h_novarCH/sqrt(2*dataEP.numberTaps);
                if scenario.flagNorm  %Normalized if set in the scenario
                    h_novarCH=h_novarCH./norm(h_novarCH);
                end
                h=h_novarCH+randn(size(h_novarCH))*sqrt(dataEP.varnoiseCH)+1i*randn(size(h_novarCH))*sqrt(dataEP.varnoiseCH);
            end
        else
            h_novarCH=dataEP.channel_h(n_h,:);
            h=h_novarCH+randn(size(h_novarCH))*sqrt(dataEP.varnoiseCH);
        end
        % Computing matrix channel
        Hmat=obth(h,dataEP.channelBlockLength+dataEP.numberTaps-1);
        Hmat=Hmat(:,1:end-dataEP.numberTaps+1);
        
        Hmat_novarCH=obth(h_novarCH,dataEP.channelBlockLength+dataEP.numberTaps-1);
        Hmat_novarCH=Hmat_novarCH(:,1:end-dataEP.numberTaps+1);

    end
    
    
    %% Initializating BCJR if scenario.flagOptimal
    if scenario.flagOptimal && ~flagMIMO
        [convmat_trans,convmat_salidas,transp_llegan_a_q,...
            transq_llegan_a_p]=var_trellis(A,dataEP.numberTaps,h.');
    end
    
    % Initialize BER for the channel
    % Hard errors
    n_errors_hard_MMSE=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_MMSE_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    
    n_mseC_MMSE_turbo=zeros(length(SNRdB),dataEP.numberTurbo);    
    n_mseD_MMSE_turbo=zeros(length(SNRdB),dataEP.numberTurbo);    
    n_mseL_MMSE_turbo=zeros(length(SNRdB),dataEP.numberTurbo);     
    
    n_errors_hard_MMSEAproxGS=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_MMSEAproxGS_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_hard_MMSEAproxNeuman=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_MMSEAproxNeuman_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_hard_MMSEAproxPCG=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_MMSEAproxPCG_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_hard_FMMSE=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_FMMSE_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_hard_BEP=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_BEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_hard_PBEP=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_PBEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);

        %%
    n_mseC_PBEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);    
    n_mseD_PBEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);    
    n_mseL_PBEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo); 
    
    n_errors_hard_DBEP=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_DBEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    
    %%
    n_mseC_DBEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);    
    n_mseD_DBEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);    
    n_mseL_DBEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);    

    
    n_errors_hard_DBEPAproxGS=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_DBEPAproxGS_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_hard_DBEPAproxNeuman=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_DBEPAproxNeuman_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_hard_DBEPAproxPCG=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_DBEPAproxPCG_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_hard_PFEP=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_PFEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_hard_SEP=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_SEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_hard_EPICLMMSE=zeros(length(SNRdB),1); % number of errors for each SNR
    n_errors_hard_EPICLMMSE_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    
    n_mseC_EPIC_turbo=zeros(length(SNRdB),dataEP.numberTurbo);    
    n_mseD_EPIC_turbo=zeros(length(SNRdB),dataEP.numberTurbo);    
    n_mseL_EPIC_turbo=zeros(length(SNRdB),dataEP.numberTurbo);    
    
    n_errors_hard_optimal=zeros(length(SNRdB),1);
    n_errors_hard_optimal_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    
    % Soft errors
    n_errors_soft_MMSEAproxGS=zeros(length(SNRdB),1);
    n_errors_soft_MMSEAproxGS_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_MMSEAproxNeuman=zeros(length(SNRdB),1);
    n_errors_soft_MMSEAproxNeuman_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_MMSEAproxPCG=zeros(length(SNRdB),1);
    n_errors_soft_MMSEAproxPCG_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_MMSE=zeros(length(SNRdB),1);
    n_errors_soft_MMSE_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_FMMSE=zeros(length(SNRdB),1);
    n_errors_soft_FMMSE_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_BEP=zeros(length(SNRdB),1);
    n_errors_soft_BEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_PBEP=zeros(length(SNRdB),1);
    n_errors_soft_PBEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_DBEP=zeros(length(SNRdB),1);
    n_errors_soft_DBEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_DBEPAproxGS=zeros(length(SNRdB),1);
    n_errors_soft_DBEPAproxGS_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_DBEPAproxNeuman=zeros(length(SNRdB),1);
    n_errors_soft_DBEPAproxNeuman_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_DBEPAproxPCG=zeros(length(SNRdB),1);
    n_errors_soft_DBEPAproxPCG_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_PFEP=zeros(length(SNRdB),1);
    n_errors_soft_PFEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_DFEP=zeros(length(SNRdB),1);
    n_errors_soft_DFEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_SEP=zeros(length(SNRdB),1);
    n_errors_soft_SEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_PKSEP=zeros(length(SNRdB),1);
    n_errors_soft_PKSEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_DKSEP=zeros(length(SNRdB),1);
    n_errors_soft_DKSEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_EPICLMMSE=zeros(length(SNRdB),1);
    n_errors_soft_EPICLMMSE_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_BPEP=zeros(length(SNRdB),1);
    n_errors_soft_BPEP_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    n_errors_soft_optimal=zeros(length(SNRdB),1);
    n_errors_soft_optimal_turbo=zeros(length(SNRdB),dataEP.numberTurbo);
    
    if dataEP.flagIterDec
        % Number of iterations of the channel decoder
        niterdec_MMSEAproxGS=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_MMSEAproxNeuman=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_MMSEAproxPCG=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_MMSE=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_FMMSE=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_BEP=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_PBEP=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_DBEP=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_DBEPAproxGS=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_DBEPAproxNeuman=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_DBEPAproxPCG=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_PFEP=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_DFEP=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_SEP=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_PKSEP=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_DKSEP=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_EPICLMMSE=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_BPEP=zeros(length(SNRdB),dataEP.numberTurbo+1);
        niterdec_Optimal=zeros(length(SNRdB),dataEP.numberTurbo+1);
    end
    
    %% Running Montecarlo for the channel and the Eb/N0s
    for snr=1:length(SNRdB)
        
        % Hard errors
        err_hard_MMSE=0;
        err_hard_MMSEturbo=zeros(1,dataEP.numberTurbo);
        err_hard_MMSEAproxGS=0;
        err_hard_MMSEAproxGSturbo=zeros(1,dataEP.numberTurbo);
        err_hard_MMSEAproxNeuman=0;
        err_hard_MMSEAproxNeumanturbo=zeros(1,dataEP.numberTurbo);
        err_hard_MMSEAproxPCG=0;
        err_hard_MMSEAproxPCGturbo=zeros(1,dataEP.numberTurbo);
        err_hard_FMMSE=0;
        err_hard_FMMSEturbo=zeros(1,dataEP.numberTurbo);
        err_hard_BEP=0;
        err_hard_BEPturbo=zeros(1,dataEP.numberTurbo);
        err_hard_PBEP=0;
        err_hard_PBEPturbo=zeros(1,dataEP.numberTurbo);
        err_hard_DBEP=0;
        err_hard_DBEPturbo=zeros(1,dataEP.numberTurbo);
        err_hard_DBEPAproxGS=0;
        err_hard_DBEPAproxGSturbo=zeros(1,dataEP.numberTurbo);
        err_hard_DBEPAproxNeuman=0;
        err_hard_DBEPAproxNeumanturbo=zeros(1,dataEP.numberTurbo);
        err_hard_DBEPAproxPCG=0;
        err_hard_DBEPAproxPCGturbo=zeros(1,dataEP.numberTurbo);
        err_hard_PFEP=0;
        err_hard_PFEPturbo=zeros(1,dataEP.numberTurbo);
        err_hard_DFEP=0;
        err_hard_DFEPturbo=zeros(1,dataEP.numberTurbo);
        err_hard_SEP=0;
        err_hard_SEPturbo=zeros(1,dataEP.numberTurbo);
        err_hard_EPICLMMSE=0;
        err_hard_EPICLMMSEturbo=zeros(1,dataEP.numberTurbo);
        err_hard_Optimal=0;
        err_hard_Optimalturbo=zeros(1,dataEP.numberTurbo);
        
        % Soft errors
        err_soft_MMSE=0;
        err_soft_MMSE_turbo=zeros(1,dataEP.numberTurbo);
        mseCit_MMSE_turbo=zeros(1,dataEP.numberTurbo);
        mseDit_MMSE_turbo=zeros(1,dataEP.numberTurbo);
        mseLit_MMSE_turbo=zeros(1,dataEP.numberTurbo);
        
        err_soft_MMSEAproxGS=0;
        err_soft_MMSEAproxGS_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_MMSEAproxNeuman=0;
        err_soft_MMSEAproxNeuman_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_MMSEAproxPCG=0;
        err_soft_MMSEAproxPCG_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_FMMSE=0;
        err_soft_FMMSE_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_BEP=0;
        err_soft_BEP_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_PBEP=0;
        err_soft_PBEP_turbo=zeros(1,dataEP.numberTurbo);
        mseCit_PBEP_turbo=zeros(1,dataEP.numberTurbo);
        mseDit_PBEP_turbo=zeros(1,dataEP.numberTurbo);
        mseLit_PBEP_turbo=zeros(1,dataEP.numberTurbo);
        
        err_soft_DBEP=0;
        err_soft_DBEP_turbo=zeros(1,dataEP.numberTurbo);
        mseCit_DBEP_turbo=zeros(1,dataEP.numberTurbo);
        mseDit_DBEP_turbo=zeros(1,dataEP.numberTurbo);
        mseLit_DBEP_turbo=zeros(1,dataEP.numberTurbo);
        
        err_soft_DBEPAproxGS=0;
        err_soft_DBEPAproxGS_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_DBEPAproxNeuman=0;
        err_soft_DBEPAproxNeuman_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_DBEPAproxPCG=0;
        err_soft_DBEPAproxPCG_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_PFEP=0;
        err_soft_PFEP_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_DFEP=0;
        err_soft_DFEP_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_SEP=0;
        err_soft_SEP_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_PKSEP=0;
        err_soft_PKSEP_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_DKSEP=0;
        err_soft_DKSEP_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_BPEP=0;
        err_soft_BPEP_turbo=zeros(1,dataEP.numberTurbo);
        err_soft_EPICLMMSE=0;
        err_soft_EPICLMMSE_turbo=zeros(1,dataEP.numberTurbo);
        mseCit_EPIC_turbo=zeros(1,dataEP.numberTurbo);
        mseDit_EPIC_turbo=zeros(1,dataEP.numberTurbo);
        mseLit_EPIC_turbo=zeros(1,dataEP.numberTurbo);
        
        err_soft_Optimal=0;
        err_soft_Optimal_turbo=zeros(1,dataEP.numberTurbo);
        
        if verbose
            disp(['n_h=',num2str(n_h),'/',num2str(dataEP.numberChannels),'  SNRdB=',num2str(SNRdB(snr)),'/',num2str(SNRdB(end))])
            %disp(['SNRdB=',num2str(SNRdB(snr)),'/',num2str(SNRdB(end))])
            %fprintf('SNRdB= %f, Clock: ',SNRdB(snr), num2str(toc(tStart)))
%            disp(['Clock: ',num2str(toc(tStart))])
            
            count=fwrite(fidLog,'n_h=','char');
            count=fwrite(fidLog,num2str(n_h),'char');
            count=fwrite(fidLog,' ','char');
            count=fwrite(fidLog,'Clock=','char');
            count=fwrite(fidLog,num2str(toc(tStart)),'char');
            count=fwrite(fidLog,10,'uint8');
        end
        
        for k=1:dataEP.numberFrames
            
            %% Sample of transmitted data
            
            %rand('state',sum(100*clock));
            %randn('state',sum(100*clock));
            %rng(simulationNumber*1e6+k,'twister');
            rng(n_h+(dataEP.numberChannels*(k-1))+(dataEP.numberChannels*dataEP.numberFrames)*(simulationNumber-1),'twister');              
            %rng(100*k*simulationNumber,'twister')
            %             rand('state',sum(100*k*n_h*dataEP.simulationNumber*snr));
            %             randn('state',sum(100*k*n_h*dataEP.simulationNumber*snr));
            
            x_sincod = randi([0 1],encoderInputLenght,1);
            
            % Channel Encoding
            x = step(hEnc, x_sincod);
            %x=[x;randi(2,[encoderOutputLenght_new-encoderOutputLenght 1])-1];
            x=[x;zeros(encoderOutputLenght_new-encoderOutputLenght,1)];
            
            
            % Interlaving
            if dataEP.flagInterleaving
                x = randintrlv(x,seed_inter);
            end
            
            %% Modulation
            x=step(hMod, x);
            
            if flagMIMO
                Limk1=modulatorOutputLengthNew/dataEP.numberAntennas(1);
                xmat=reshape(x,dataEP.numberAntennas(1),Limk1);
            else
                Limk1=modulatorOutputLengthNew/dataEP.channelBlockLength;
                xmat=reshape(x,dataEP.channelBlockLength,Limk1);
            end
            
            % We keep the vector of noise and received vector to use the
            % same ones when turbo is computed
            noiseALL=zeros(size(Hmat,1),Limk1);
            rALLC=zeros(size(Hmat,1),Limk1);
            
            % For double proposals (it is needed if the frame is divided
            % into different blocks)
            meanE_DBEP=zeros(size(xmat));
            varE_DBEP=zeros(size(xmat));
            meanE_DBEPAproxGS=zeros(size(xmat));
            varE_DBEPAproxGS=zeros(size(xmat));
            meanE_DBEPAproxNeuman=zeros(size(xmat));
            varE_DBEPAproxNeuman=zeros(size(xmat));
            meanE_DBEPAproxPCG=zeros(size(xmat));
            varE_DBEPAproxPCG=zeros(size(xmat));
            meanE_DFEP=zeros(size(xmat));
            varE_DFEP=zeros(size(xmat));
            meanE_DKSEP=zeros(size(xmat));
            varE_DKSEP=zeros(size(xmat));
            meanE_BPEP=zeros(size(xmat));
            varE_BPEP=zeros(size(xmat));
            
            for idxTurbo=0:dataEP.numberTurbo
                LLR_MMSE=[];
                LLR_MMSEAproxGS=[];
                LLR_MMSEAproxNeuman=[];
                LLR_MMSEAproxPCG=[];
                LLR_FMMSE=[];
                LLR_BEP=[];
                LLR_EPICLMMSE=[];
                LLR_SEP=[];
                LLR_PBEP=[];
                LLR_DBEP=[];
                LLR_DBEPAproxGS=[];
                LLR_DBEPAproxNeuman=[];
                LLR_DBEPAproxPCG=[];
                LLR_PFEP=[];
                LLR_DFEP=[];
                LLR_PKSEP=[];
                LLR_DKSEP=[];
                LLR_BPEP=[];
                LLR_Optimal=[];
                
                for k1=1:Limk1
                    x=xmat(:,k1); %each column of xmat is a vector transmitted by the Tx antennas
                    
                    xinitial=x;
                    if idxTurbo==0
                        % AWGN noise
                        noise=sigma_novarCH(snr)*randn(size(Hmat_novarCH,1),1);
                        if dataEP.complexFlag
                            noise=noise+1i*sigma_novarCH(snr)*randn(size(Hmat_novarCH,1),1);
                        end
                        noiseALL(:,k1)=noise;
                        
                        % Salida del canal
                        Hx=Hmat_novarCH*xinitial; %JJMF
                        rcomplex=Hmat_novarCH*xinitial+noise;
                        rALLC(:,k1)=rcomplex;
                    else
                        noise=noiseALL(:,k1);
                        rcomplex=rALLC(:,k1);
                    end
                    
            %%%%%%%%%%%%%
            %% DETECTION
            %%%%%%%%%%%%%%
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%% MMSE [Muranov10] %%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagMMSE
                        if scenario.SNRdBMMSE(1)<=SNRdB(snr) && scenario.SNRdBMMSE(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_MMSE,prob_b_MMSE]=MMSEalg(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[]);
                                % Hard detection
                                err_hard_MMSE=err_hard_MMSE+sum(xinitial~=x_decod_hard_MMSE);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_MMSE,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_MMSE=[LLR_MMSE real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_MMSE,prob_b_MMSE]=MMSEalg(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_ldpc_MMSE(:,:,k1).');
                                % Hard detection
                                err_hard_MMSEturbo(idxTurbo)=err_hard_MMSEturbo(idxTurbo)+sum(xinitial~=x_decod_hard_MMSE);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_MMSE,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_MMSE=[LLR_MMSE real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%% approximated MMSE with GS %%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagMMSEAproxGS
                        if scenario.SNRdBMMSE(1)<=SNRdB(snr) && scenario.SNRdBMMSE(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_MMSEAproxGS,prob_b_MMSEAproxGS]=MMSEalgAproxGS(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[]);
                                % Hard detection
                                err_hard_MMSEAproxGS=err_hard_MMSEAproxGS+sum(xinitial~=x_decod_hard_MMSEAproxGS);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_MMSEAproxGS,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_MMSEAproxGS=[LLR_MMSEAproxGS real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_MMSEAproxGS,prob_b_MMSEAproxGS]=MMSEalgAproxGS(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_ldpc_MMSEAproxGS(:,:,k1).');
                                % Hard detection
                                err_hard_MMSEAproxGSturbo(idxTurbo)=err_hard_MMSEAproxGSturbo(idxTurbo)+sum(xinitial~=x_decod_hard_MMSEAproxGS);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_MMSEAproxGS,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_MMSEAproxGS=[LLR_MMSEAproxGS real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%% approximated MMSE with Neuman %%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagMMSEAproxNeuman
                        if scenario.SNRdBMMSE(1)<=SNRdB(snr) && scenario.SNRdBMMSE(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_MMSEAproxNeuman,prob_b_MMSEAproxNeuman]=MMSEalgAproxNeuman(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[]);
                                % Hard detection
                                err_hard_MMSEAproxNeuman=err_hard_MMSEAproxNeuman+sum(xinitial~=x_decod_hard_MMSEAproxNeuman);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_MMSEAproxNeuman,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_MMSEAproxNeuman=[LLR_MMSEAproxNeuman real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_MMSEAproxNeuman,prob_b_MMSEAproxNeuman]=MMSEalgAproxNeuman(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_ldpc_MMSEAproxNeuman(:,:,k1).');
                                % Hard detection
                                err_hard_MMSEAproxNeumanturbo(idxTurbo)=err_hard_MMSEAproxNeumanturbo(idxTurbo)+sum(xinitial~=x_decod_hard_MMSEAproxNeuman);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_MMSEAproxNeuman,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_MMSEAproxNeuman=[LLR_MMSEAproxNeuman real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%% approximated MMSE with PCG %%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagMMSEAproxPCG
                        if scenario.SNRdBMMSE(1)<=SNRdB(snr) && scenario.SNRdBMMSE(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_MMSEAproxPCG,prob_b_MMSEAproxPCG]=MMSEalgAproxPCG(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[]);
                                % Hard detection
                                err_hard_MMSEAproxPCG=err_hard_MMSEAproxPCG+sum(xinitial~=x_decod_hard_MMSEAproxPCG);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_MMSEAproxPCG,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_MMSEAproxPCG=[LLR_MMSEAproxPCG real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_MMSEAproxPCG,prob_b_MMSEAproxPCG]=MMSEalgAproxPCG(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_ldpc_MMSEAproxPCG(:,:,k1).');
                                % Hard detection
                                err_hard_MMSEAproxPCGturbo(idxTurbo)=err_hard_MMSEAproxPCGturbo(idxTurbo)+sum(xinitial~=x_decod_hard_MMSEAproxPCG);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_MMSEAproxPCG,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_MMSEAproxPCG=[LLR_MMSEAproxPCG real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%% F-MMSE (Wiener) [Tuchler02] %%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagFMMSE
                        if scenario.SNRdBMMSE(1)<=SNRdB(snr) && scenario.SNRdBMMSE(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_FMMSE,prob_b_FMMSE]=FMMSEalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,Hmat,[]);
                                % Hard detection
                                err_hard_FMMSE=err_hard_FMMSE+sum(xinitial~=x_decod_hard_FMMSE);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_FMMSE,symbolmap,dataEP.flagInterleaving,seed_inter);
                                %prob_b_after_demap=prob_b_BEP;
                                LLR_FMMSE=[LLR_FMMSE real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_FMMSE,prob_b_FMMSE]=FMMSEalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,Hmat,pui_ldpc_FMMSE(:,:,k1).');
                                % Hard detection
                                err_hard_FMMSEturbo(idxTurbo)=err_hard_FMMSEturbo(idxTurbo)+sum(xinitial~=x_decod_hard_FMMSE);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_FMMSE,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_FMMSE=[LLR_FMMSE real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% BEP [Santos16] %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagBEP
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_BEP,prob_b_BEP]=BEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[],dataEP.BEP_S);
                                % Hard detection
                                err_hard_BEP=err_hard_BEP+sum(xinitial~=x_decod_hard_BEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_BEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_BEP=[LLR_BEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_BEP,prob_b_BEP]=BEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_ldpc_BEP(:,:,k1).',dataEP.BEP_S);
                                % Hard detection
                                err_hard_BEPturbo(idxTurbo)=err_hard_BEPturbo(idxTurbo)+sum(xinitial~=x_decod_hard_BEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_BEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_BEP=[LLR_BEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% P-BEP [Santos18] %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagPBEP
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_PBEP,prob_b_PBEP]=PBEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[],idxTurbo,dataEP.PBEP_S);
                                % Hard detection
                                err_hard_PBEP=err_hard_PBEP+sum(xinitial~=x_decod_hard_PBEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_PBEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_PBEP=[LLR_PBEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_PBEP,prob_b_PBEP,nu_qi,nu_piDec,nu_piDecEP]=PBEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_ldpc_PBEP(:,:,k1).',idxTurbo,dataEP.PBEP_S);
                                % Hard detection
                                err_hard_PBEPturbo(idxTurbo)=err_hard_PBEPturbo(idxTurbo)+sum(xinitial~=x_decod_hard_PBEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_PBEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_PBEP=[LLR_PBEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                                    mseCit_PBEP_turbo(idxTurbo)=sum(abs(nu_piDec-xinitial).^2)/energy/length(Hx);
                                    mseLit_PBEP_turbo(idxTurbo)=sum(abs(nu_piDecEP-xinitial).^2)/energy/length(Hx);                                   
                                    mseDit_PBEP_turbo(idxTurbo)=sum(abs(nu_qi-xinitial).^2)/energy/length(Hx);                               
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% D-BEP [Santos19]  %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagDBEP
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_DBEP,prob_b_DBEP,meanE_DBEPk1,varE_DBEPk1]=DBEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[],idxTurbo,[],[],dataEP.DBEP_S);
                                meanE_DBEP(:,k1)=meanE_DBEPk1;
                                varE_DBEP(:,k1)=varE_DBEPk1;
                                % Hard detection
                                err_hard_DBEP=err_hard_DBEP+sum(xinitial~=x_decod_hard_DBEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DBEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DBEP=[LLR_DBEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                                %[idxTurbo, SNRdB(snr), sum(abs(Hmat*meanE_DBEPk1-Hx).^2)/energy/length(Hx)]
                            else
                                [x_decod_hard_DBEP,prob_b_DBEP,meanE_DBEPk1,varE_DBEPk1,nu_qi,nu_piDec,nu_piDecEP]=DBEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_ldpc_DBEP(:,:,k1).',idxTurbo,meanE_DBEP(:,k1),varE_DBEP(:,k1),dataEP.DBEP_S);
                                meanE_DBEP(:,k1)=meanE_DBEPk1;
                                varE_DBEP(:,k1)=varE_DBEPk1;
                                % Hard detection
                                err_hard_DBEPturbo(idxTurbo)=err_hard_DBEPturbo(idxTurbo)+sum(xinitial~=x_decod_hard_DBEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DBEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DBEP=[LLR_DBEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
%                                 [idxTurbo, SNRdB(snr), sum(abs(Hmat*nu_qi-Hx).^2)/energy/length(Hx), ...
%                                     ...%sum(abs(Hmat*meanE_DBEPk1-Hx).^2)/energy/length(Hx)...
%                                     sum(abs(nu_piDec-xinitial).^2)/energy/length(Hx)...
%                                     sum(abs(nu_piDecEP-xinitial).^2)/energy/length(Hx)...                                    
%                                     sum(abs(nu_qi-xinitial).^2)/energy/length(Hx)]
                               
                                    mseCit_DBEP_turbo(idxTurbo)=sum(abs(nu_piDec-xinitial).^2)/energy/length(Hx);
                                    mseLit_DBEP_turbo(idxTurbo)=sum(abs(nu_piDecEP-xinitial).^2)/energy/length(Hx);                                   
                                    mseDit_DBEP_turbo(idxTurbo)=sum(abs(nu_qi-xinitial).^2)/energy/length(Hx);
                                
%                                Hqi=Hmat*nu_qi;
%                                [Hqi(1:10,:) Hx(1:10,:)]
%                                 energy
%                                [xinitial(1:10),nu_piDec(1:10), nu_piDecEP(1:10), nu_qi(1:10)]
                            end
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% approximated D-BEP with GS %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagDBEPAproxGS
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_DBEPAproxGS,prob_b_DBEPAproxGS,meanE_DBEPAproxk1,varE_DBEPAproxk1]=DBEPalgAproxGS(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[],idxTurbo,[],[]);
                                meanE_DBEPAproxGS(:,k1)=meanE_DBEPAproxk1;
                                varE_DBEPAproxGS(:,k1)=varE_DBEPAproxk1;
                                % Hard detection
                                err_hard_DBEPAproxGS=err_hard_DBEPAproxGS+sum(xinitial~=x_decod_hard_DBEPAproxGS);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DBEPAproxGS,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DBEPAproxGS=[LLR_DBEPAproxGS real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_DBEPAproxGS,prob_b_DBEPAproxGS,meanE_DBEPAproxk1,varE_DBEPAproxk1]=DBEPalgAproxGS(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_ldpc_DBEPAproxGS(:,:,k1).',idxTurbo,meanE_DBEPAproxGS(:,k1),varE_DBEPAproxGS(:,k1));
                                meanE_DBEPAproxGS(:,k1)=meanE_DBEPAproxk1;
                                varE_DBEPAproxGS(:,k1)=varE_DBEPAproxk1;
                                % Hard detection
                                err_hard_DBEPAproxGSturbo(idxTurbo)=err_hard_DBEPAproxGSturbo(idxTurbo)+sum(xinitial~=x_decod_hard_DBEPAproxGS);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DBEPAproxGS,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DBEPAproxGS=[LLR_DBEPAproxGS real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% approximated D-BEP with Neuman %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagDBEPAproxNeuman
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_DBEPAproxNeuman,prob_b_DBEPAproxNeuman,meanE_DBEPAproxk1,varE_DBEPAproxk1]=DBEPalgAproxNeuman(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[],idxTurbo,[],[]);
                                meanE_DBEPAproxNeuman(:,k1)=meanE_DBEPAproxk1;
                                varE_DBEPAproxNeuman(:,k1)=varE_DBEPAproxk1;
                                % Hard detection
                                err_hard_DBEPAproxNeuman=err_hard_DBEPAproxNeuman+sum(xinitial~=x_decod_hard_DBEPAproxNeuman);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DBEPAproxNeuman,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DBEPAproxNeuman=[LLR_DBEPAproxNeuman real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_DBEPAproxNeuman,prob_b_DBEPAproxNeuman,meanE_DBEPAproxk1,varE_DBEPAproxk1]=DBEPalgAproxNeuman(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_ldpc_DBEPAproxNeuman(:,:,k1).',idxTurbo,meanE_DBEPAproxNeuman(:,k1),varE_DBEPAproxNeuman(:,k1));
                                meanE_DBEPAproxNeuman(:,k1)=meanE_DBEPAproxk1;
                                varE_DBEPAproxNeuman(:,k1)=varE_DBEPAproxk1;
                                % Hard detection
                                err_hard_DBEPAproxNeumanturbo(idxTurbo)=err_hard_DBEPAproxNeumanturbo(idxTurbo)+sum(xinitial~=x_decod_hard_DBEPAproxNeuman);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DBEPAproxNeuman,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DBEPAproxNeuman=[LLR_DBEPAproxNeuman real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% approximated D-BEP with PCG %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagDBEPAproxPCG
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_DBEPAproxPCG,prob_b_DBEPAproxPCG,meanE_DBEPAproxk1,varE_DBEPAproxk1]=DBEPalgAproxPCG(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[],idxTurbo,[],[]);
                                meanE_DBEPAproxPCG(:,k1)=meanE_DBEPAproxk1;
                                varE_DBEPAproxPCG(:,k1)=varE_DBEPAproxk1;
                                % Hard detection
                                err_hard_DBEPAproxPCG=err_hard_DBEPAproxPCG+sum(xinitial~=x_decod_hard_DBEPAproxPCG);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DBEPAproxPCG,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DBEPAproxPCG=[LLR_DBEPAproxPCG real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_DBEPAproxPCG,prob_b_DBEPAproxPCG,meanE_DBEPAproxk1,varE_DBEPAproxk1]=DBEPalgAproxPCG(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_ldpc_DBEPAproxPCG(:,:,k1).',idxTurbo,meanE_DBEPAproxPCG(:,k1),varE_DBEPAproxPCG(:,k1));
                                meanE_DBEPAproxPCG(:,k1)=meanE_DBEPAproxk1;
                                varE_DBEPAproxPCG(:,k1)=varE_DBEPAproxk1;
                                % Hard detection
                                err_hard_DBEPAproxPCGturbo(idxTurbo)=err_hard_DBEPAproxPCGturbo(idxTurbo)+sum(xinitial~=x_decod_hard_DBEPAproxPCG);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DBEPAproxPCG,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DBEPAproxPCG=[LLR_DBEPAproxPCG real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% P-FEP [Santos18] %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagPFEP
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_PFEP,prob_b_PFEP]=PFEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,Hmat,[],idxTurbo);
                                % Hard detection
                                err_hard_PFEP=err_hard_PFEP+sum(xinitial~=x_decod_hard_PFEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_PFEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_PFEP=[LLR_PFEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_PFEP,prob_b_PFEP]=PFEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,Hmat,pui_ldpc_PFEP(:,:,k1).',idxTurbo);
                                % Hard detection
                                err_hard_PFEPturbo(idxTurbo)=err_hard_PFEPturbo(idxTurbo)+sum(xinitial~=x_decod_hard_PFEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_PFEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_PFEP=[LLR_PFEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% D-FEP [Santos19] %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagDFEP
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_DFEP,prob_b_DFEP,~,varE_DFEPk1,meanE_DFEPk1]=DFEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,Hmat,[],idxTurbo,[],[]);
                                meanE_DFEP(:,k1)=meanE_DFEPk1;
                                varE_DFEP(:,k1)=varE_DFEPk1;
                                % Hard detection
                                err_hard_DFEP=err_hard_DFEP+sum(xinitial~=x_decod_hard_DFEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DFEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DFEP=[LLR_DFEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_DFEP,prob_b_DFEP,~,varE_DFEPk1,meanE_DFEPk1]=DFEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,Hmat,pui_ldpc_DFEP(:,:,k1).',idxTurbo,meanE_DFEP(:,k1),varE_DFEP(:,k1));
                                meanE_DFEP(:,k1)=meanE_DFEPk1;
                                varE_DFEP(:,k1)=varE_DFEPk1;
                                % Hard detection
                                err_hard_DFEPturbo(idxTurbo)=err_hard_DFEPturbo(idxTurbo)+sum(xinitial~=x_decod_hard_DFEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DFEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DFEP=[LLR_DFEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% SEP [Santos16] %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagSEP
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_SEP,prob_b_SEP]=SEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,[]);
                                % Hard detection
                                err_hard_SEP=err_hard_SEP+sum(xinitial~=x_decod_hard_SEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_SEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_SEP=[LLR_SEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_SEP,prob_b_SEP]=SEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,pui_ldpc_SEP(:,:,k1).');
                                % Hard detection
                                err_hard_SEPturbo(idxTurbo)=err_hard_SEPturbo(idxTurbo)+sum(xinitial~=x_decod_hard_SEP);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_SEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_SEP=[LLR_SEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% P-KSEP [Santos18c] %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagPKSEP
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                %[~,prob_b_PKSEP]=PKSEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,[],idxTurbo);
                                [~,prob_b_PKSEP]=PKSEPalg_OptComp(A,dataEP.complexFlag,sigma(snr),rcomplex,h,[],idxTurbo,dataEP.KSEP_S);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_PKSEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_PKSEP=[LLR_PKSEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                %[~,prob_b_PKSEP]=PKSEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,pui_ldpc_PKSEP(:,:,k1).',idxTurbo);
                                [~,prob_b_PKSEP]=PKSEPalg_OptComp(A,dataEP.complexFlag,sigma(snr),rcomplex,h,pui_ldpc_PKSEP(:,:,k1).',idxTurbo,dataEP.KSEP_S);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_PKSEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_PKSEP=[LLR_PKSEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% D-KSEP [Santos19] %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagDKSEP
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                %[~,prob_b_DKSEP,meanE_DKSEPk1,varE_DKSEPk1]=DKSEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,[],idxTurbo,[],[]);
                                [~,prob_b_DKSEP,meanE_DKSEPk1,varE_DKSEPk1]=DKSEPalg_OptComp(A,dataEP.complexFlag,sigma(snr),rcomplex,h,[],idxTurbo,[],[],dataEP.DKSEP_S);
%                                k1
%                                [size(meanE_DKSEPk1)' size(meanE_DKSEP)']
                                meanE_DKSEP(:,k1)=meanE_DKSEPk1;
                                varE_DKSEP(:,k1)=varE_DKSEPk1;
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DKSEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DKSEP=[LLR_DKSEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                %[~,prob_b_DKSEP,meanE_DKSEPk1,varE_DKSEPk1]=DKSEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,pui_ldpc_DKSEP(:,:,k1).',idxTurbo,meanE_DKSEP(:,k1),varE_DKSEP(:,k1));
                                [~,prob_b_DKSEP,meanE_DKSEPk1,varE_DKSEPk1]=DKSEPalg_OptComp(A,dataEP.complexFlag,sigma(snr),rcomplex,h,pui_ldpc_DKSEP(:,:,k1).',idxTurbo,meanE_DKSEP(:,k1),varE_DKSEP(:,k1),dataEP.DKSEP_S);
                                meanE_DKSEP(:,k1)=meanE_DKSEPk1;
                                varE_DKSEP(:,k1)=varE_DKSEPk1;
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_DKSEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_DKSEP=[LLR_DKSEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% BP-EP [Sun15] %%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagBPEP
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [~,prob_b_BPEP,meanE_BPEPk1,varE_BPEPk1]=BPEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,[],[],[]);
                                meanE_BPEP(:,k1)=meanE_BPEPk1;
                                varE_BPEP(:,k1)=varE_BPEPk1;
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_BPEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_BPEP=[LLR_BPEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [~,prob_b_BPEP,meanE_BPEPk1,varE_BPEPk1]=BPEPalg(A,dataEP.complexFlag,sigma(snr),rcomplex,h,pui_ldpc_BPEP(:,:,k1).',meanE_BPEP(:,k1),varE_BPEP(:,k1));
                                meanE_BPEP(:,k1)=meanE_BPEPk1;
                                varE_BPEP(:,k1)=varE_BPEPk1;
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_BPEP,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_BPEP=[LLR_BPEP real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% EP IC-LMMSE [Senst11] %%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagEPICLMMSE
                        if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                            if idxTurbo==0
                                [x_decod_hard_EPICLMMSE,prob_b_EPICLMMSE]=EPICLMMSEalg(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[]);
                                % Hard detection
                                err_hard_EPICLMMSE=err_hard_EPICLMMSE+sum(xinitial~=x_decod_hard_EPICLMMSE);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_EPICLMMSE,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_EPICLMMSE=[LLR_EPICLMMSE real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            else
                                [x_decod_hard_EPICLMMSE,prob_b_EPICLMMSE]=EPICLMMSEalg(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_ldpc_EPICLMMSE(:,:,k1).');
                                % Hard detection
                                err_hard_EPICLMMSEturbo(idxTurbo)=err_hard_EPICLMMSEturbo(idxTurbo)+sum(xinitial~=x_decod_hard_EPICLMMSE);
                                % Demodulation
                                prob_b_after_demap=demap(prob_b_EPICLMMSE,symbolmap,dataEP.flagInterleaving,seed_inter);
                                LLR_EPICLMMSE=[LLR_EPICLMMSE real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                            end
                        end
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%% BCJR or MAP %%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if scenario.flagOptimal
                        if scenario.SNRdBOptimal(1)<=SNRdB(snr) && scenario.SNRdBOptimal(2)>=SNRdB(snr)
                            if flagMIMO %MIMO: ML
                                if idxTurbo==0
                                    [x_decod_hard_Optimal,prob_b_Optimal]=ML(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,[]);
                                    % Hard detection
                                    err_hard_Optimal=err_hard_Optimal+sum(xinitial~=x_decod_hard_Optimal);
                                    % Demodulation
                                    prob_b_after_demap=demap(prob_b_Optimal,symbolmap,dataEP.flagInterleaving,seed_inter);
                                    LLR_Optimal=[LLR_Optimal real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                                else
                                    [x_decod_hard_Optimal,prob_b_Optimal]=ML(A,dataEP.complexFlag,sigma(snr),rcomplex,Hmat,pui_Optimal(:,:,k1));
                                    % Hard detection
                                    err_hard_Optimalturbo(idxTurbo)=err_hard_Optimalturbo(idxTurbo)+sum(xinitial~=x_decod_hard_Optimal);
                                    % Demodulation
                                    prob_b_after_demap=demap(prob_b_Optimal,symbolmap,dataEP.flagInterleaving,seed_inter);
                                    LLR_Optimal=[LLR_Optimal real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                                    
                                    
                                end
                            else %EQUALIZATION: BCJR
                                if idxTurbo==0
                                    [x_decod_hard_Optimal,prob_b_Optimal]=BCJRalg(A,dataEP.complexFlag,dataEP.channelBlockLength+dataEP.numberTaps-1,h,[xinitial; A(1)*ones(dataEP.numberTaps-1,1)],noise,sigma(snr),convmat_trans,convmat_salidas,transp_llegan_a_q,transq_llegan_a_p,[]);
                                    % Hard detection
                                    err_hard_Optimal=err_hard_Optimal+sum(xinitial~=x_decod_hard_Optimal(1:dataEP.channelBlockLength));
                                    % Demodulation
                                    prob_b_after_demap=demap(prob_b_Optimal(:,1:dataEP.channelBlockLength),symbolmap,dataEP.flagInterleaving,seed_inter);
                                    LLR_Optimal=[LLR_Optimal real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                                else
                                    [x_decod_hard_Optimal,prob_b_Optimal]=BCJRalg(A,dataEP.complexFlag,dataEP.channelBlockLength+dataEP.numberTaps-1,h,[xinitial; A(1)*ones(dataEP.numberTaps-1,1)],noise,sigma(snr),convmat_trans,convmat_salidas,transp_llegan_a_q,transq_llegan_a_p,pui_Optimal(:,:,k1));
                                    % Hard detection
                                    err_hard_Optimalturbo(idxTurbo)=err_hard_Optimalturbo(idxTurbo)+sum(xinitial~=x_decod_hard_Optimal(1:dataEP.channelBlockLength));
                                    % Demodulation
                                    prob_b_after_demap=demap(prob_b_Optimal(:,1:dataEP.channelBlockLength),symbolmap,dataEP.flagInterleaving,seed_inter);
                                    LLR_Optimal=[LLR_Optimal real(log(prob_b_after_demap(1,:)./prob_b_after_demap(2,:)))];
                                end
                            end
                        end
                    end
               end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Computing the soft errors & DECODING
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% MMSE [Muranov10] %%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagMMSE
                    if scenario.SNRdBMMSE(1)<=SNRdB(snr) && scenario.SNRdBMMSE(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_MMSE=max(min(LLR_MMSE(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_MMSE,niterdec] = step(hDec, LLR_input_decoder_MMSE');
                        if dataEP.flagIterDec, niterdec_MMSE(snr,idxTurbo+1)=niterdec_MMSE(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_MMSE_uncoded=LLR_output_decoder_MMSE(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_MMSE=LLR_output_decoder_MMSE-LLR_input_decoder_MMSE';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_MMSE)); % p(ct=0|y)
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0; % p(ct=1|y)
                            p01=[p0';p1'];
                            pui_ldpc_MMSE=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter); % p(xt=Ai|y)
                            pui_ldpc_MMSE=reshape(pui_ldpc_MMSE,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_MMSE_uncoded)); % p(ct=0|y)
                        p1=1-p0; % p(ct=1|y)
                        % Soft detection (of the information word)
                        x_decod_soft_MMSE=round(p1);
                        if idxTurbo==0
                            err_soft_MMSE=err_soft_MMSE+sum(xor(x_sincod,x_decod_soft_MMSE));
                        else
                            err_soft_MMSE_turbo(idxTurbo)=err_soft_MMSE_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_MMSE));
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% approximated MMSE with GS %%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagMMSEAproxGS
                    if scenario.SNRdBMMSE(1)<=SNRdB(snr) && scenario.SNRdBMMSE(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_MMSEAproxGS=max(min(LLR_MMSEAproxGS(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_MMSEAproxGS,niterdec] = step(hDec, LLR_input_decoder_MMSEAproxGS');
                        if dataEP.flagIterDec, niterdec_MMSEAproxGS(snr,idxTurbo+1)=niterdec_MMSEAproxGS(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_MMSEAproxGS_uncoded=LLR_output_decoder_MMSEAproxGS(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_MMSEAproxGS=LLR_output_decoder_MMSEAproxGS-LLR_input_decoder_MMSEAproxGS';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_MMSEAproxGS)); % p(ct=0|y)
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0; % p(ct=1|y)
                            p01=[p0';p1'];
                            pui_ldpc_MMSEAproxGS=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter); % p(xt=Ai|y)
                            pui_ldpc_MMSEAproxGS=reshape(pui_ldpc_MMSEAproxGS,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_MMSEAproxGS_uncoded)); % p(ct=0|y)
                        p1=1-p0; % p(ct=1|y)
                        % Soft detection (of the information word)
                        x_decod_soft_MMSEAproxGS=round(p1);
                        if idxTurbo==0
                            err_soft_MMSEAproxGS=err_soft_MMSEAproxGS+sum(xor(x_sincod,x_decod_soft_MMSEAproxGS));
                        else
                            err_soft_MMSEAproxGS_turbo(idxTurbo)=err_soft_MMSEAproxGS_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_MMSEAproxGS));
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% approximated MMSE with Neuman %%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagMMSEAproxNeuman
                    if scenario.SNRdBMMSE(1)<=SNRdB(snr) && scenario.SNRdBMMSE(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_MMSEAproxNeuman=max(min(LLR_MMSEAproxNeuman(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_MMSEAproxNeuman,niterdec] = step(hDec, LLR_input_decoder_MMSEAproxNeuman');
                        if dataEP.flagIterDec, niterdec_MMSEAproxNeuman(snr,idxTurbo+1)=niterdec_MMSEAproxNeuman(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_MMSEAproxNeuman_uncoded=LLR_output_decoder_MMSEAproxNeuman(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_MMSEAproxNeuman=LLR_output_decoder_MMSEAproxNeuman-LLR_input_decoder_MMSEAproxNeuman';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_MMSEAproxNeuman)); % p(ct=0|y)
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0; % p(ct=1|y)
                            p01=[p0';p1'];
                            pui_ldpc_MMSEAproxNeuman=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter); % p(xt=Ai|y)
                            pui_ldpc_MMSEAproxNeuman=reshape(pui_ldpc_MMSEAproxNeuman,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_MMSEAproxNeuman_uncoded)); % p(ct=0|y)
                        p1=1-p0; % p(ct=1|y)
                        % Soft detection (of the information word)
                        x_decod_soft_MMSEAproxNeuman=round(p1);
                        if idxTurbo==0
                            err_soft_MMSEAproxNeuman=err_soft_MMSEAproxNeuman+sum(xor(x_sincod,x_decod_soft_MMSEAproxNeuman));
                        else
                            err_soft_MMSEAproxNeuman_turbo(idxTurbo)=err_soft_MMSEAproxNeuman_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_MMSEAproxNeuman));
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% approximated MMSE with PCG %%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagMMSEAproxPCG
                    if scenario.SNRdBMMSE(1)<=SNRdB(snr) && scenario.SNRdBMMSE(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_MMSEAproxPCG=max(min(LLR_MMSEAproxPCG(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_MMSEAproxPCG,niterdec] = step(hDec, LLR_input_decoder_MMSEAproxPCG');
                        if dataEP.flagIterDec, niterdec_MMSEAproxPCG(snr,idxTurbo+1)=niterdec_MMSEAproxPCG(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_MMSEAproxPCG_uncoded=LLR_output_decoder_MMSEAproxPCG(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_MMSEAproxPCG=LLR_output_decoder_MMSEAproxPCG-LLR_input_decoder_MMSEAproxPCG';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_MMSEAproxPCG)); % p(ct=0|y)
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0; % p(ct=1|y)
                            p01=[p0';p1'];
                            pui_ldpc_MMSEAproxPCG=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter); % p(xt=Ai|y)
                            pui_ldpc_MMSEAproxPCG=reshape(pui_ldpc_MMSEAproxPCG,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_MMSEAproxPCG_uncoded)); % p(ct=0|y)
                        p1=1-p0; % p(ct=1|y)
                        % Soft detection (of the information word)
                        x_decod_soft_MMSEAproxPCG=round(p1);
                        if idxTurbo==0
                            err_soft_MMSEAproxPCG=err_soft_MMSEAproxPCG+sum(xor(x_sincod,x_decod_soft_MMSEAproxPCG));
                        else
                            err_soft_MMSEAproxPCG_turbo(idxTurbo)=err_soft_MMSEAproxPCG_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_MMSEAproxPCG));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% F-MMSE (Wiener) [Tuchler02] %%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagFMMSE
                    if scenario.SNRdBMMSE(1)<=SNRdB(snr) && scenario.SNRdBMMSE(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_FMMSE=max(min(LLR_FMMSE(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_FMMSE,niterdec] = step(hDec, LLR_input_decoder_FMMSE');
                        if dataEP.flagIterDec, niterdec_FMMSE(snr,idxTurbo+1)=niterdec_FMMSE(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_FMMSE_uncoded=LLR_output_decoder_FMMSE(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_FMMSE=LLR_output_decoder_FMMSE-LLR_input_decoder_FMMSE';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_FMMSE)); % p(ct=0|y)
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0; % p(ct=1|y)
                            p01=[p0';p1'];
                            pui_ldpc_FMMSE=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter); % p(xt=Ai|y)
                            pui_ldpc_FMMSE=reshape(pui_ldpc_FMMSE,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_FMMSE_uncoded)); % p(ct=0|y)
                        p1=1-p0; % p(ct=1|y)
                        % Soft detection (of the information word)
                        x_decod_soft_FMMSE=round(p1);
                        if idxTurbo==0
                            err_soft_FMMSE=err_soft_FMMSE+sum(xor(x_sincod,x_decod_soft_FMMSE));
                        else
                            err_soft_FMMSE_turbo(idxTurbo)=err_soft_FMMSE_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_FMMSE));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% BEP [Santos16] %%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagBEP
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_BEP=max(min(LLR_BEP(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_BEP,niterdec] = step(hDec, LLR_input_decoder_BEP');
                        if dataEP.flagIterDec, niterdec_BEP(snr,idxTurbo+1)=niterdec_BEP(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_BEP_uncoded=LLR_output_decoder_BEP(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_BEP=LLR_output_decoder_BEP-LLR_input_decoder_BEP';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_BEP));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc_BEP=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_BEP=reshape(pui_ldpc_BEP,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_BEP_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_BEP=round(p1);
                        if idxTurbo==0
                            err_soft_BEP=err_soft_BEP+sum(xor(x_sincod,x_decod_soft_BEP));
                        else
                            err_soft_BEP_turbo(idxTurbo)=err_soft_BEP_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_BEP));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% P-BEP [Santos18] %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagPBEP
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_PBEP=max(min(LLR_PBEP(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_PBEP,niterdec] = step(hDec, LLR_input_decoder_PBEP');
                        if dataEP.flagIterDec, niterdec_PBEP(snr,idxTurbo+1)=niterdec_PBEP(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_PBEP_uncoded=LLR_output_decoder_PBEP(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_PBEP=LLR_output_decoder_PBEP-LLR_input_decoder_PBEP';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_PBEP));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc_PBEP=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_PBEP=reshape(pui_ldpc_PBEP,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_PBEP_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_PBEP=round(p1);
                        if idxTurbo==0
                            err_soft_PBEP=err_soft_PBEP+sum(xor(x_sincod,x_decod_soft_PBEP));
                        else
                            err_soft_PBEP_turbo(idxTurbo)=err_soft_PBEP_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_PBEP));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% D-BEP [Santos19]  %%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagDBEP
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_DBEP=max(min(LLR_DBEP(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_DBEP,niterdec] = step(hDec, LLR_input_decoder_DBEP');
                        if dataEP.flagIterDec, niterdec_DBEP(snr,idxTurbo+1)=niterdec_DBEP(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_DBEP_uncoded=LLR_output_decoder_DBEP(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_DBEP=LLR_output_decoder_DBEP-LLR_input_decoder_DBEP';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_DBEP));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc_DBEP=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_DBEP=reshape(pui_ldpc_DBEP,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_DBEP_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_DBEP=round(p1);
                        if idxTurbo==0
                            err_soft_DBEP=err_soft_DBEP+sum(xor(x_sincod,x_decod_soft_DBEP));
                        else
                            err_soft_DBEP_turbo(idxTurbo)=err_soft_DBEP_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_DBEP));
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% approximated D-BEP with GS %%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagDBEPAproxGS
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_DBEPAproxGS=max(min(LLR_DBEPAproxGS(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_DBEPAproxGS,niterdec] = step(hDec, LLR_input_decoder_DBEPAproxGS');
                        if dataEP.flagIterDec, niterdec_DBEPAproxGS(snr,idxTurbo+1)=niterdec_DBEPAproxGS(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_DBEPAproxGS_uncoded=LLR_output_decoder_DBEPAproxGS(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_DBEPAproxGS=LLR_output_decoder_DBEPAproxGS-LLR_input_decoder_DBEPAproxGS';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_DBEPAproxGS));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc_DBEPAproxGS=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_DBEPAproxGS=reshape(pui_ldpc_DBEPAproxGS,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_DBEPAproxGS_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_DBEPAproxGS=round(p1);
                        if idxTurbo==0
                            err_soft_DBEPAproxGS=err_soft_DBEPAproxGS+sum(xor(x_sincod,x_decod_soft_DBEPAproxGS));
                        else
                            err_soft_DBEPAproxGS_turbo(idxTurbo)=err_soft_DBEPAproxGS_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_DBEPAproxGS));
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% approximated D-BEP with Neuman %%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagDBEPAproxNeuman
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_DBEPAproxNeuman=max(min(LLR_DBEPAproxNeuman(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_DBEPAproxNeuman,niterdec] = step(hDec, LLR_input_decoder_DBEPAproxNeuman');
                        if dataEP.flagIterDec, niterdec_DBEPAproxNeuman(snr,idxTurbo+1)=niterdec_DBEPAproxNeuman(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_DBEPAproxNeuman_uncoded=LLR_output_decoder_DBEPAproxNeuman(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_DBEPAproxNeuman=LLR_output_decoder_DBEPAproxNeuman-LLR_input_decoder_DBEPAproxNeuman';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_DBEPAproxNeuman));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc_DBEPAproxNeuman=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_DBEPAproxNeuman=reshape(pui_ldpc_DBEPAproxNeuman,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_DBEPAproxNeuman_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_DBEPAproxNeuman=round(p1);
                        if idxTurbo==0
                            err_soft_DBEPAproxNeuman=err_soft_DBEPAproxNeuman+sum(xor(x_sincod,x_decod_soft_DBEPAproxNeuman));
                        else
                            err_soft_DBEPAproxNeuman_turbo(idxTurbo)=err_soft_DBEPAproxNeuman_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_DBEPAproxNeuman));
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% approximated D-BEP with PCG %%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagDBEPAproxPCG
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_DBEPAproxPCG=max(min(LLR_DBEPAproxPCG(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_DBEPAproxPCG,niterdec] = step(hDec, LLR_input_decoder_DBEPAproxPCG');
                        if dataEP.flagIterDec, niterdec_DBEPAproxPCG(snr,idxTurbo+1)=niterdec_DBEPAproxPCG(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_DBEPAproxPCG_uncoded=LLR_output_decoder_DBEPAproxPCG(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_DBEPAproxPCG=LLR_output_decoder_DBEPAproxPCG-LLR_input_decoder_DBEPAproxPCG';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_DBEPAproxPCG));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc_DBEPAproxPCG=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_DBEPAproxPCG=reshape(pui_ldpc_DBEPAproxPCG,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_DBEPAproxPCG_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_DBEPAproxPCG=round(p1);
                        if idxTurbo==0
                            err_soft_DBEPAproxPCG=err_soft_DBEPAproxPCG+sum(xor(x_sincod,x_decod_soft_DBEPAproxPCG));
                        else
                            err_soft_DBEPAproxPCG_turbo(idxTurbo)=err_soft_DBEPAproxPCG_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_DBEPAproxPCG));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% P-FEP [Santos18] %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagPFEP
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_PFEP=max(min(LLR_PFEP(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_PFEP,niterdec] = step(hDec, LLR_input_decoder_PFEP');
                        if dataEP.flagIterDec, niterdec_PFEP(snr,idxTurbo+1)=niterdec_PFEP(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_PFEP_uncoded=LLR_output_decoder_PFEP(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_PFEP=LLR_output_decoder_PFEP-LLR_input_decoder_PFEP';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_PFEP));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc_PFEP=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_PFEP=reshape(pui_ldpc_PFEP,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_PFEP_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_PFEP=round(p1);
                        if idxTurbo==0
                            err_soft_PFEP=err_soft_PFEP+sum(xor(x_sincod,x_decod_soft_PFEP));
                        else
                            err_soft_PFEP_turbo(idxTurbo)=err_soft_PFEP_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_PFEP));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% D-FEP [Santos19] %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagDFEP
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_DFEP=max(min(LLR_DFEP(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_DFEP,niterdec] = step(hDec, LLR_input_decoder_DFEP');
                        if dataEP.flagIterDec, niterdec_DFEP(snr,idxTurbo+1)=niterdec_DFEP(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_DFEP_uncoded=LLR_output_decoder_DFEP(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_DFEP=LLR_output_decoder_DFEP-LLR_input_decoder_DFEP';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_DFEP));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc_DFEP=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_DFEP=reshape(pui_ldpc_DFEP,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_DFEP_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_DFEP=round(p1);
                        if idxTurbo==0
                            err_soft_DFEP=err_soft_DFEP+sum(xor(x_sincod,x_decod_soft_DFEP));
                        else
                            err_soft_DFEP_turbo(idxTurbo)=err_soft_DFEP_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_DFEP));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% SEP [Santos16] %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagSEP
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_SEP=max(min(LLR_SEP(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_SEP,niterdec] = step(hDec, LLR_input_decoder_SEP');
                        if dataEP.flagIterDec, niterdec_SEP(snr,idxTurbo+1)=niterdec_SEP(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_SEP_uncoded=LLR_output_decoder_SEP(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_SEP=LLR_output_decoder_SEP-LLR_input_decoder_SEP';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_SEP));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc_SEP=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_SEP=reshape(pui_ldpc_SEP,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_SEP_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_SEP=round(p1);
                        if idxTurbo==0
                            err_soft_SEP=err_soft_SEP+sum(xor(x_sincod,x_decod_soft_SEP));
                        else
                            err_soft_SEP_turbo(idxTurbo)=err_soft_SEP_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_SEP));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% P-KSEP [Santos18c] %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagPKSEP
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_PKSEP=max(min(LLR_PKSEP(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_PKSEP,niterdec] = step(hDec, LLR_input_decoder_PKSEP');
                        if dataEP.flagIterDec, niterdec_PKSEP(snr,idxTurbo+1)=niterdec_PKSEP(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_PKSEP_uncoded=LLR_output_decoder_PKSEP(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_PKSEP=LLR_output_decoder_PKSEP-LLR_input_decoder_PKSEP';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_PKSEP));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_PKSEP=reshape(pui_ldpc,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_PKSEP_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_PKSEP=round(p1);
                        if idxTurbo==0
                            err_soft_PKSEP=err_soft_PKSEP+sum(xor(x_sincod,x_decod_soft_PKSEP));
                        else
                            err_soft_PKSEP_turbo(idxTurbo)=err_soft_PKSEP_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_PKSEP));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% D-KSEP [Santos19] %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagDKSEP
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_DKSEP=max(min(LLR_DKSEP(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_DKSEP,niterdec] = step(hDec, LLR_input_decoder_DKSEP');
                        if dataEP.flagIterDec, niterdec_DKSEP(snr,idxTurbo+1)=niterdec_DKSEP(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_DKSEP_uncoded=LLR_output_decoder_DKSEP(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_DKSEP=LLR_output_decoder_DKSEP-LLR_input_decoder_DKSEP';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_DKSEP));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_DKSEP=reshape(pui_ldpc,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_DKSEP_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_DKSEP=round(p1);
                        if idxTurbo==0
                            err_soft_DKSEP=err_soft_DKSEP+sum(xor(x_sincod,x_decod_soft_DKSEP));
                        else
                            err_soft_DKSEP_turbo(idxTurbo)=err_soft_DKSEP_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_DKSEP));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% BP-EP [Sun15] %%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagBPEP
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_BPEP=max(min(LLR_BPEP(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_BPEP,niterdec] = step(hDec, LLR_input_decoder_BPEP');
                        if dataEP.flagIterDec, niterdec_BPEP(snr,idxTurbo+1)=niterdec_BPEP(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_BPEP_uncoded=LLR_output_decoder_BPEP(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_BPEP=LLR_output_decoder_BPEP-LLR_input_decoder_BPEP';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_BPEP));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_BPEP=reshape(pui_ldpc,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_BPEP_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_BPEP=round(p1);
                        if idxTurbo==0
                            err_soft_BPEP=err_soft_BPEP+sum(xor(x_sincod,x_decod_soft_BPEP));
                        else
                            err_soft_BPEP_turbo(idxTurbo)=err_soft_BPEP_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_BPEP));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% EP IC-LMMSE [Senst11] %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagEPICLMMSE
                    if scenario.SNRdBEP(1)<=SNRdB(snr) && scenario.SNRdBEP(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr
                        LLR_input_decoder_EPICLMMSE=max(min(LLR_EPICLMMSE(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_EPICLMMSE,niterdec] = step(hDec, LLR_input_decoder_EPICLMMSE');
                        if dataEP.flagIterDec, niterdec_EPICLMMSE(snr,idxTurbo+1)=niterdec_EPICLMMSE(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_EPICLMMSE_uncoded=LLR_output_decoder_EPICLMMSE(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_EPICLMMSE=LLR_output_decoder_EPICLMMSE-LLR_input_decoder_EPICLMMSE';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_EPICLMMSE));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            pui_ldpc_EPICLMMSE=reshape(pui_ldpc,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_EPICLMMSE_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_EPICLMMSE=round(p1);
                        if idxTurbo==0
                            err_soft_EPICLMMSE=err_soft_EPICLMMSE+sum(xor(x_sincod,x_decod_soft_EPICLMMSE));
                        else
                            err_soft_EPICLMMSE_turbo(idxTurbo)=err_soft_EPICLMMSE_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_EPICLMMSE));
                        end
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% BCJR or MAP %%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if scenario.flagOptimal
                    if scenario.SNRdBOptimal(1)<=SNRdB(snr) && scenario.SNRdBOptimal(2)>=SNRdB(snr)
                        % We limit the maximum and minimum value of llr to
                        % avoid instabilities
                        LLR_input_decoder_Optimal=max(min(LLR_Optimal(:,1:encoderOutputLenght),limLDPC),-limLDPC);
                        
                        % Posterior at the output of the LDPC Decoder
                        [LLR_output_decoder_Optimal,niterdec] = step(hDec, LLR_input_decoder_Optimal');
                        if dataEP.flagIterDec, niterdec_Optimal(snr,idxTurbo+1)=niterdec_Optimal(snr,idxTurbo+1)+niterdec; end
                        % Posterior uncoded
                        LLR_output_decoder_Optimal_uncoded=LLR_output_decoder_Optimal(1:encoderInputLenght);
                        % Extrinsic (at the output of the LDPC Decoder)
                        LLR_output_decoder_Optimal=LLR_output_decoder_Optimal-LLR_input_decoder_Optimal';
                        
                        if dataEP.numberTurbo>0
                            p0=1./(1+exp(-LLR_output_decoder_Optimal));
                            p0=[p0;ones(encoderOutputLenght_new-encoderOutputLenght,1)];
                            p1=1-p0;
                            p01=[p0';p1'];
                            pui_ldpc=map(p01,symbolmap,dataEP.flagInterleaving,seed_inter);
                            
                            pui_Optimal=pui_ldpc;
                            pui_Optimal=reshape(pui_Optimal,dataEP.M,modulatorOutputLengthNew/Limk1,Limk1);
                            %pui_Optimal=pui_Optimal-LLR_Optimal;
                        end
                        
                        p0=1./(1+exp(-LLR_output_decoder_Optimal_uncoded));
                        p1=1-p0;
                        % Soft detection
                        x_decod_soft_Optimal=round(p1);
                        if idxTurbo==0
                            err_soft_Optimal=err_soft_Optimal+sum(xor(x_sincod,x_decod_soft_Optimal));
                        else
                            err_soft_Optimal_turbo(idxTurbo)=err_soft_Optimal_turbo(idxTurbo)+sum(xor(x_sincod,x_decod_soft_Optimal));
                        end
                    end
                end
            end
        end
        
        % Hard errors
        n_errors_hard_MMSE(snr)=n_errors_hard_MMSE(snr)+err_hard_MMSE;
        n_errors_hard_MMSE_turbo(snr,:)=n_errors_hard_MMSE_turbo(snr,:)+err_hard_MMSEturbo;
        n_errors_hard_MMSEAproxGS(snr)=n_errors_hard_MMSEAproxGS(snr)+err_hard_MMSEAproxGS;
        n_errors_hard_MMSEAproxGS_turbo(snr,:)=n_errors_hard_MMSEAproxGS_turbo(snr,:)+err_hard_MMSEAproxGSturbo;
        n_errors_hard_MMSEAproxNeuman(snr)=n_errors_hard_MMSEAproxNeuman(snr)+err_hard_MMSEAproxNeuman;
        n_errors_hard_MMSEAproxNeuman_turbo(snr,:)=n_errors_hard_MMSEAproxNeuman_turbo(snr,:)+err_hard_MMSEAproxNeumanturbo;
        n_errors_hard_MMSEAproxPCG(snr)=n_errors_hard_MMSEAproxPCG(snr)+err_hard_MMSEAproxPCG;
        n_errors_hard_MMSEAproxPCG_turbo(snr,:)=n_errors_hard_MMSEAproxPCG_turbo(snr,:)+err_hard_MMSEAproxPCGturbo;
        n_errors_hard_FMMSE(snr)=n_errors_hard_FMMSE(snr)+err_hard_FMMSE;
        n_errors_hard_FMMSE_turbo(snr,:)=n_errors_hard_FMMSE_turbo(snr,:)+err_hard_FMMSEturbo;
        n_errors_hard_BEP(snr)=n_errors_hard_BEP(snr)+err_hard_BEP;
        n_errors_hard_BEP_turbo(snr,:)=n_errors_hard_BEP_turbo(snr,:)+err_hard_BEPturbo;
        n_errors_hard_PBEP(snr)=n_errors_hard_PBEP(snr)+err_hard_PBEP;
        n_errors_hard_PBEP_turbo(snr,:)=n_errors_hard_PBEP_turbo(snr,:)+err_hard_PBEPturbo;
        n_errors_hard_DBEP(snr)=n_errors_hard_DBEP(snr)+err_hard_DBEP;
        n_errors_hard_DBEP_turbo(snr,:)=n_errors_hard_DBEP_turbo(snr,:)+err_hard_DBEPturbo;
        n_errors_hard_DBEPAproxGS(snr)=n_errors_hard_DBEPAproxGS(snr)+err_hard_DBEPAproxGS;
        n_errors_hard_DBEPAproxGS_turbo(snr,:)=n_errors_hard_DBEPAproxGS_turbo(snr,:)+err_hard_DBEPAproxGSturbo;
        n_errors_hard_DBEPAproxNeuman(snr)=n_errors_hard_DBEPAproxNeuman(snr)+err_hard_DBEPAproxNeuman;
        n_errors_hard_DBEPAproxNeuman_turbo(snr,:)=n_errors_hard_DBEPAproxNeuman_turbo(snr,:)+err_hard_DBEPAproxNeumanturbo;
        n_errors_hard_DBEPAproxPCG(snr)=n_errors_hard_DBEPAproxPCG(snr)+err_hard_DBEPAproxPCG;
        n_errors_hard_DBEPAproxPCG_turbo(snr,:)=n_errors_hard_DBEPAproxPCG_turbo(snr,:)+err_hard_DBEPAproxPCGturbo;
        n_errors_hard_PFEP(snr)=n_errors_hard_PFEP(snr)+err_hard_PFEP;
        n_errors_hard_PFEP_turbo(snr,:)=n_errors_hard_PFEP_turbo(snr,:)+err_hard_PFEPturbo;
        n_errors_hard_SEP(snr)=n_errors_hard_SEP(snr)+err_hard_SEP;
        n_errors_hard_SEP_turbo(snr,:)=n_errors_hard_SEP_turbo(snr,:)+err_hard_SEPturbo;
        n_errors_hard_EPICLMMSE(snr)=n_errors_hard_EPICLMMSE(snr)+err_hard_EPICLMMSE;
        n_errors_hard_EPICLMMSE_turbo(snr,:)=n_errors_hard_EPICLMMSE_turbo(snr,:)+err_hard_EPICLMMSEturbo;
        n_errors_hard_optimal(snr)=n_errors_hard_optimal(snr)+err_hard_Optimal;
        n_errors_hard_optimal_turbo(snr,:)=n_errors_hard_optimal_turbo(snr,:)+err_hard_Optimalturbo;
        
        % Soft errors
        n_errors_soft_MMSE(snr)=n_errors_soft_MMSE(snr)+err_soft_MMSE;
        n_errors_soft_MMSE_turbo(snr,:)=n_errors_soft_MMSE_turbo(snr,:)+err_soft_MMSE_turbo;
        n_errors_soft_MMSEAproxGS(snr)=n_errors_soft_MMSEAproxGS(snr)+err_soft_MMSEAproxGS;
        n_errors_soft_MMSEAproxGS_turbo(snr,:)=n_errors_soft_MMSEAproxGS_turbo(snr,:)+err_soft_MMSEAproxGS_turbo;
        n_errors_soft_MMSEAproxNeuman(snr)=n_errors_soft_MMSEAproxNeuman(snr)+err_soft_MMSEAproxNeuman;
        n_errors_soft_MMSEAproxNeuman_turbo(snr,:)=n_errors_soft_MMSEAproxNeuman_turbo(snr,:)+err_soft_MMSEAproxNeuman_turbo;
        n_errors_soft_MMSEAproxPCG(snr)=n_errors_soft_MMSEAproxPCG(snr)+err_soft_MMSEAproxPCG;
        n_errors_soft_MMSEAproxPCG_turbo(snr,:)=n_errors_soft_MMSEAproxPCG_turbo(snr,:)+err_soft_MMSEAproxPCG_turbo;
        n_errors_soft_FMMSE(snr)=n_errors_soft_FMMSE(snr)+err_soft_FMMSE;
        n_errors_soft_FMMSE_turbo(snr,:)=n_errors_soft_FMMSE_turbo(snr,:)+err_soft_FMMSE_turbo;
        n_errors_soft_BEP(snr)=n_errors_soft_BEP(snr)+err_soft_BEP;
        n_errors_soft_BEP_turbo(snr,:)=n_errors_soft_BEP_turbo(snr,:)+err_soft_BEP_turbo;
        n_errors_soft_PBEP(snr)=n_errors_soft_PBEP(snr)+err_soft_PBEP;
        n_errors_soft_PBEP_turbo(snr,:)=n_errors_soft_PBEP_turbo(snr,:)+err_soft_PBEP_turbo;

        n_mseC_PBEP_turbo(snr,:)=n_mseC_PBEP_turbo(snr,:)+mseCit_PBEP_turbo;
        n_mseD_PBEP_turbo(snr,:)=n_mseD_PBEP_turbo(snr,:)+mseDit_PBEP_turbo;
        n_mseL_PBEP_turbo(snr,:)=n_mseL_PBEP_turbo(snr,:)+mseLit_PBEP_turbo;
        
        n_errors_soft_DBEP(snr)=n_errors_soft_DBEP(snr)+err_soft_DBEP;
        n_errors_soft_DBEP_turbo(snr,:)=n_errors_soft_DBEP_turbo(snr,:)+err_soft_DBEP_turbo;
        
        n_mseC_DBEP_turbo(snr,:)=n_mseC_DBEP_turbo(snr,:)+mseCit_DBEP_turbo;
        n_mseD_DBEP_turbo(snr,:)=n_mseD_DBEP_turbo(snr,:)+mseDit_DBEP_turbo;
        n_mseL_DBEP_turbo(snr,:)=n_mseL_DBEP_turbo(snr,:)+mseLit_DBEP_turbo;
        
        n_errors_soft_DBEPAproxGS(snr)=n_errors_soft_DBEPAproxGS(snr)+err_soft_DBEPAproxGS;
        n_errors_soft_DBEPAproxGS_turbo(snr,:)=n_errors_soft_DBEPAproxGS_turbo(snr,:)+err_soft_DBEPAproxGS_turbo;
        n_errors_soft_DBEPAproxNeuman(snr)=n_errors_soft_DBEPAproxNeuman(snr)+err_soft_DBEPAproxNeuman;
        n_errors_soft_DBEPAproxNeuman_turbo(snr,:)=n_errors_soft_DBEPAproxNeuman_turbo(snr,:)+err_soft_DBEPAproxNeuman_turbo;
        n_errors_soft_DBEPAproxPCG(snr)=n_errors_soft_DBEPAproxPCG(snr)+err_soft_DBEPAproxPCG;
        n_errors_soft_DBEPAproxPCG_turbo(snr,:)=n_errors_soft_DBEPAproxPCG_turbo(snr,:)+err_soft_DBEPAproxPCG_turbo;
        n_errors_soft_PFEP(snr)=n_errors_soft_PFEP(snr)+err_soft_PFEP;
        n_errors_soft_PFEP_turbo(snr,:)=n_errors_soft_PFEP_turbo(snr,:)+err_soft_PFEP_turbo;
        n_errors_soft_DFEP(snr)=n_errors_soft_DFEP(snr)+err_soft_DFEP;
        n_errors_soft_DFEP_turbo(snr,:)=n_errors_soft_DFEP_turbo(snr,:)+err_soft_DFEP_turbo;
        n_errors_soft_SEP(snr)=n_errors_soft_SEP(snr)+err_soft_SEP;
        n_errors_soft_SEP_turbo(snr,:)=n_errors_soft_SEP_turbo(snr,:)+err_soft_SEP_turbo;
        n_errors_soft_PKSEP(snr)=n_errors_soft_PKSEP(snr)+err_soft_PKSEP;
        n_errors_soft_PKSEP_turbo(snr,:)=n_errors_soft_PKSEP_turbo(snr,:)+err_soft_PKSEP_turbo;
        n_errors_soft_DKSEP(snr)=n_errors_soft_DKSEP(snr)+err_soft_DKSEP;
        n_errors_soft_DKSEP_turbo(snr,:)=n_errors_soft_DKSEP_turbo(snr,:)+err_soft_DKSEP_turbo;
        n_errors_soft_BPEP(snr)=n_errors_soft_BPEP(snr)+err_soft_BPEP;
        n_errors_soft_BPEP_turbo(snr,:)=n_errors_soft_BPEP_turbo(snr,:)+err_soft_BPEP_turbo;
        n_errors_soft_EPICLMMSE(snr)=n_errors_soft_EPICLMMSE(snr)+err_soft_EPICLMMSE;
        n_errors_soft_EPICLMMSE_turbo(snr,:)=n_errors_soft_EPICLMMSE_turbo(snr,:)+err_soft_EPICLMMSE_turbo;
        n_errors_soft_optimal(snr)=n_errors_soft_optimal(snr)+err_soft_Optimal;
        n_errors_soft_optimal_turbo(snr,:)=n_errors_soft_optimal_turbo(snr,:)+err_soft_Optimal_turbo;
    end
    
    % Hard BER
    BER_hard_MMSE=n_errors_hard_MMSE/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_MMSE_turbo=n_errors_hard_MMSE_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_MMSEAproxGS=n_errors_hard_MMSEAproxGS/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_MMSEAproxGS_turbo=n_errors_hard_MMSEAproxGS_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_MMSEAproxNeuman=n_errors_hard_MMSEAproxNeuman/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_MMSEAproxNeuman_turbo=n_errors_hard_MMSEAproxNeuman_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_MMSEAproxPCG=n_errors_hard_MMSEAproxPCG/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_MMSEAproxPCG_turbo=n_errors_hard_MMSEAproxPCG_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_FMMSE=n_errors_hard_FMMSE/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_FMMSE_turbo=n_errors_hard_FMMSE_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_BEP=n_errors_hard_BEP/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_BEP_turbo=n_errors_hard_BEP_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_PBEP=n_errors_hard_PBEP/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_PBEP_turbo=n_errors_hard_PBEP_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_DBEP=n_errors_hard_DBEP/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_DBEP_turbo=n_errors_hard_DBEP_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_DBEPAproxGS=n_errors_hard_DBEPAproxGS/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_DBEPAproxGS_turbo=n_errors_hard_DBEPAproxGS_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_DBEPAproxNeuman=n_errors_hard_DBEPAproxNeuman/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_DBEPAproxNeuman_turbo=n_errors_hard_DBEPAproxNeuman_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_DBEPAproxPCG=n_errors_hard_DBEPAproxPCG/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_DBEPAproxPCG_turbo=n_errors_hard_DBEPAproxPCG_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_PFEP=n_errors_hard_PFEP/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_PFEP_turbo=n_errors_hard_PFEP_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_SEP=n_errors_hard_SEP/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_SEP_turbo=n_errors_hard_SEP_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_EPICLMMSE=n_errors_hard_EPICLMMSE/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_EPICLMMSE_turbo=n_errors_hard_EPICLMMSE_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_Optimal=n_errors_hard_optimal/(modulatorOutputLengthNew*dataEP.numberFrames);
    BER_hard_Optimal_turbo=n_errors_hard_optimal_turbo/(modulatorOutputLengthNew*dataEP.numberFrames);
    
    % Soft BER
    BER_soft_MMSE=n_errors_soft_MMSE/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_MMSE_turbo=n_errors_soft_MMSE_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_MMSEAproxGS=n_errors_soft_MMSEAproxGS/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_MMSEAproxGS_turbo=n_errors_soft_MMSEAproxGS_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_MMSEAproxNeuman=n_errors_soft_MMSEAproxNeuman/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_MMSEAproxNeuman_turbo=n_errors_soft_MMSEAproxNeuman_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_MMSEAproxPCG=n_errors_soft_MMSEAproxPCG/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_MMSEAproxPCG_turbo=n_errors_soft_MMSEAproxPCG_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_FMMSE=n_errors_soft_FMMSE/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_FMMSE_turbo=n_errors_soft_FMMSE_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_BEP=n_errors_soft_BEP/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_BEP_turbo=n_errors_soft_BEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_PBEP=n_errors_soft_PBEP/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_PBEP_turbo=n_errors_soft_PBEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    
    MSED_PBEP_turbo=n_mseD_PBEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    MSEC_PBEP_turbo=n_mseC_PBEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    MSEL_PBEP_turbo=n_mseL_PBEP_turbo/(encoderInputLenght*dataEP.numberFrames);
       
    BER_soft_DBEP=n_errors_soft_DBEP/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_DBEP_turbo=n_errors_soft_DBEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    
    MSED_DBEP_turbo=n_mseD_DBEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    MSEC_DBEP_turbo=n_mseC_DBEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    MSEL_DBEP_turbo=n_mseL_DBEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    
    BER_soft_DBEPAproxGS=n_errors_soft_DBEPAproxGS/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_DBEPAproxGS_turbo=n_errors_soft_DBEPAproxGS_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_DBEPAproxNeuman=n_errors_soft_DBEPAproxNeuman/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_DBEPAproxNeuman_turbo=n_errors_soft_DBEPAproxNeuman_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_DBEPAproxPCG=n_errors_soft_DBEPAproxPCG/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_DBEPAproxPCG_turbo=n_errors_soft_DBEPAproxPCG_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_PFEP=n_errors_soft_PFEP/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_PFEP_turbo=n_errors_soft_PFEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_DFEP=n_errors_soft_DFEP/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_DFEP_turbo=n_errors_soft_DFEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_SEP=n_errors_soft_SEP/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_SEP_turbo=n_errors_soft_SEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_PKSEP=n_errors_soft_PKSEP/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_PKSEP_turbo=n_errors_soft_PKSEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_DKSEP=n_errors_soft_DKSEP/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_DKSEP_turbo=n_errors_soft_DKSEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_BPEP=n_errors_soft_BPEP/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_BPEP_turbo=n_errors_soft_BPEP_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_EPICLMMSE=n_errors_soft_EPICLMMSE/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_EPICLMMSE_turbo=n_errors_soft_EPICLMMSE_turbo/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_Optimal=n_errors_soft_optimal/(encoderInputLenght*dataEP.numberFrames);
    BER_soft_Optimal_turbo=n_errors_soft_optimal_turbo/(encoderInputLenght*dataEP.numberFrames);
    
    if dataEP.flagIterDec
        % Number of iterations of the channel decoder
        niterdec_MMSEAproxGS=niterdec_MMSEAproxGS/dataEP.numberFrames;
        niterdec_MMSEAproxNeuman=niterdec_MMSEAproxNeuman/dataEP.numberFrames;
        niterdec_MMSEAproxPCG=niterdec_MMSEAproxPCG/dataEP.numberFrames;
        niterdec_MMSE=niterdec_MMSE/dataEP.numberFrames;
        niterdec_FMMSE=niterdec_FMMSE/dataEP.numberFrames;
        niterdec_BEP=niterdec_BEP/dataEP.numberFrames;
        niterdec_PBEP=niterdec_PBEP/dataEP.numberFrames;
        niterdec_DBEP=niterdec_DBEP/dataEP.numberFrames;
        niterdec_DBEPAproxGS=niterdec_DBEPAproxGS/dataEP.numberFrames;
        niterdec_DBEPAproxNeuman=niterdec_DBEPAproxNeuman/dataEP.numberFrames;
        niterdec_DBEPAproxPCG=niterdec_DBEPAproxPCG/dataEP.numberFrames;
        niterdec_PFEP=niterdec_PFEP/dataEP.numberFrames;
        niterdec_DFEP=niterdec_DFEP/dataEP.numberFrames;
        niterdec_SEP=niterdec_SEP/dataEP.numberFrames;
        niterdec_PKSEP=niterdec_PKSEP/dataEP.numberFrames;
        niterdec_DKSEP=niterdec_DKSEP/dataEP.numberFrames;
        niterdec_EPICLMMSE=niterdec_EPICLMMSE/dataEP.numberFrames;
        niterdec_BPEP=niterdec_BPEP/dataEP.numberFrames;
        niterdec_Optimal=niterdec_Optimal/dataEP.numberFrames;
    end
    
    fidaux=fopen(nameO,'a');
    fidauxMSE=fopen(nameMSE,'a');
    count=fwrite(fidaux,num2str(n_h),'char');
    count=fwrite(fidaux,' ','char');
    count=fwrite(fidauxMSE,num2str(n_h),'char');
    count=fwrite(fidauxMSE,' ','char');
    for snr=1:length(SNRdB)
        count=fwrite(fidaux,num2str(BER_hard_MMSE(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_MMSE_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_MMSEAproxGS(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_MMSEAproxGS_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_MMSEAproxNeuman(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_MMSEAproxNeuman_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_MMSEAproxPCG(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_MMSEAproxPCG_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_FMMSE(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_FMMSE_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_BEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_BEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_PBEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_PBEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_DBEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_DBEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_DBEPAproxGS(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_DBEPAproxGS_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_DBEPAproxNeuman(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_DBEPAproxNeuman_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_DBEPAproxPCG(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_DBEPAproxPCG_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_PFEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_PFEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_SEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_SEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_EPICLMMSE(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_EPICLMMSE_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_hard_Optimal(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_hard_Optimal_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        
        count=fwrite(fidaux,num2str(BER_soft_MMSE(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_MMSE_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_MMSEAproxGS(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_MMSEAproxGS_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_MMSEAproxNeuman(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_MMSEAproxNeuman_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_MMSEAproxPCG(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_MMSEAproxPCG_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_FMMSE(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_FMMSE_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_BEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_BEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_PBEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_PBEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_DBEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_DBEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_DBEPAproxGS(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_DBEPAproxGS_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_DBEPAproxNeuman(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_DBEPAproxNeuman_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_DBEPAproxPCG(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_DBEPAproxPCG_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_PFEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_PFEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_DFEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_DFEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_SEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_SEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_PKSEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_PKSEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_DKSEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_DKSEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_BPEP(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_BPEP_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_EPICLMMSE(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_EPICLMMSE_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        count=fwrite(fidaux,num2str(BER_soft_Optimal(snr)),'char');
        count=fwrite(fidaux,' ','uint8');
        for k1=1:dataEP.numberTurbo
            count=fwrite(fidaux,num2str(BER_soft_Optimal_turbo(snr,k1)),'char');
            count=fwrite(fidaux,' ','uint8');
        end
        if dataEP.flagIterDec
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_MMSEAproxGS(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_MMSEAproxNeuman(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_MMSEAproxPCG(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_MMSE(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_FMMSE(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_BEP(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_PBEP(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_DBEP(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_DBEPAproxGS(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_DBEPAproxNeuman(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_DBEPAproxPCG(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_PFEP(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_DFEP(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_SEP(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_PKSEP(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_DKSEP(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_EPICLMMSE(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_BPEP(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            for k1=1:dataEP.numberTurbo+1
                count=fwrite(fidaux,num2str(niterdec_Optimal(snr,k1)),'char');
                count=fwrite(fidaux,' ','uint8');
            end
            
      
                  
        end
        %%MSE 
            %DBEP
            %count=fwrite(fidauxMSE,num2str(MSEC_DBEP_turbo),'char');
            %count=fwrite(fidauxMSE,' ','uint8');        
            for k1=1:dataEP.numberTurbo
                count=fwrite(fidauxMSE,num2str(MSEC_DBEP_turbo(snr,k1)),'char');
                count=fwrite(fidauxMSE,' ','uint8');
            end
            %count=fwrite(fidauxMSE,num2str(MSED_DBEP_turbo),'char');
            %count=fwrite(fidauxMSE,' ','uint8');
            for k1=1:dataEP.numberTurbo
                count=fwrite(fidauxMSE,num2str(MSED_DBEP_turbo(snr,k1)),'char');
                count=fwrite(fidauxMSE,' ','uint8');
            end
            %count=fwrite(fidauxMSE,num2str(MSEL_DBEP_turbo),'char');
            %count=fwrite(fidauxMSE,' ','uint8');
            for k1=1:dataEP.numberTurbo
                count=fwrite(fidauxMSE,num2str(MSEL_DBEP_turbo(snr,k1)),'char');
                count=fwrite(fidauxMSE,' ','uint8');
            end      
            %count=fwrite(fidauxMSE,num2str(MSEC_DBEP_turbo),'char');
            %count=fwrite(fidauxMSE,' ','uint8');  
            %%PBEP
            for k1=1:dataEP.numberTurbo
                count=fwrite(fidauxMSE,num2str(MSEC_PBEP_turbo(snr,k1)),'char');
                count=fwrite(fidauxMSE,' ','uint8');
            end
            %count=fwrite(fidauxMSE,num2str(MSED_DBEP_turbo),'char');
            %count=fwrite(fidauxMSE,' ','uint8');
            for k1=1:dataEP.numberTurbo
                count=fwrite(fidauxMSE,num2str(MSED_PBEP_turbo(snr,k1)),'char');
                count=fwrite(fidauxMSE,' ','uint8');
            end
            %count=fwrite(fidauxMSE,num2str(MSEL_DBEP_turbo),'char');
            %count=fwrite(fidauxMSE,' ','uint8');
            for k1=1:dataEP.numberTurbo
                count=fwrite(fidauxMSE,num2str(MSEL_PBEP_turbo(snr,k1)),'char');
                count=fwrite(fidauxMSE,' ','uint8');
            end    
    end
    fwrite(fidaux,10,'uint8');
    fwrite(fidauxMSE,10,'uint8');
    %Note that every row of the output file is the result ofr a channel.
    %To average them sum the rows and divide by numberChannels
    fclose(fidaux);
    fclose(fidauxMSE);
    
end

%% Closing files
%fclose(fidaux);
if verbose
    count=fwrite(fidLog,'End','char');
    fclose(fidLog);
end

