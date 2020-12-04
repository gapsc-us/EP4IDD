% Select Equalization or MIMO
dataEP.scenario=2; % 1 for 'Equalization' or 2 for 'MIMO'

% Monte Carlo is used, with numberSimulations experiments with a
% transmission of numberFrames coded words
dataEP.numberSimulations=20; %number of diferent simulations, or average using MC
dataEP.numberFrames=1;  %number of words transmitted per simulation

% Limit to the LLRs
dataEP.LLRlim=3;

% Modulation
dataEP.M=256; %Modulation order
dataEP.complexFlag=1;   %If 1 the constellation is complex
dataEP.flagPSK=0; % if complexFlag==0, this flag will be ignored. If 1, the
% complex constellation is a PSK, while if 0 we will use a QAM


% EbN0 range
% dataEP.SNRdBIni=5-3; %11.5-3;
 dataEP.SNRdBStep=.5;%.5;
% dataEP.SNRdBEnd=20-3; %15-3
%Different EbN0 can be selected for some methods
scenario.SNRdBEP=[dataEP.SNRdBIni,dataEP.SNRdBEnd];
scenario.SNRdBMMSE=[dataEP.SNRdBIni,dataEP.SNRdBEnd]; 
scenario.SNRdBOptimal=[dataEP.SNRdBIni,dataEP.SNRdBEnd]; 
 

%EP inner iterations
    dataEP.BEP_S=10;% BEP [Santos16], Default 10
    dataEP.PBEP_S=3; % P-BEP [Santos18], Default 3
    dataEP.DBEP_S=1; % D-BEP [Santos19], Default 1
    dataEP.KSEP_S=3;    
    dataEP.DKSEP_S=1; % D-KSEP [Santos19], Default 1
    dataEP.flagEPiterations=0; %If 1, same value will be used for all. 
    %This value must be set in dataEP.EPiterations when calling mainEP
    %approaches
    
 %outer iterations (turbo = IDD) 
 dataEP.numberTurbo=5; % number of turbo iterations. If no turbo (standalone
                      %equalization), set it to 0.

% LDCP channel coded block length
    dataEP.long=1;
    %dataEP.rate=1/4;
    dataEP.channelBlockLength=64800;
    %dataEP.bits=dataEP.channelBlockLength
    dataEP.bits=floor(log2(dataEP.channelBlockLength*dataEP.rate)); %rate 1/2
 
    %dataEP.channelBlockLength=512;
    %dataEP.bits=log2(dataEP.channelBlockLength)-1; %rate 1/2
    dataEP.flagBlockLengths=1; %Must be set to 1 if analyzing over Block Length
    dataEP.maxIterLDPC=20; %maximun number of BP iterations of the LDPC decoder
       % the decoder stops if all parities are satisfied
    
    %n=<bits. If n< 2^bits/rate the transmitted block is splitted into
    %(2^bits/rate)/n transmissions, each one with initiation and termination
    dataEP.flagInterleaving=0; % To 1 if interleaving and deinterleaving is included.
    %Set it to 0 in other case. Note that in LDCP this has no impact.


%% Defining the figures to be plotted
dataEP.flagOnlySoft=1; % To 1 if we only plot Soft estimations
dataEP.flagIterDec=0; % To 1 if we also plot the number of iterations of the channel decoder
dataEP.flagBERalongT=0; % To 1 if we also plot the BER versus the number of turbo iterations


% Channel paramenters
dataEP.channel_h=[]; %[ 1 2 3 2 1]/sqrt(19); % channel_h=[] means random channels. We can also set a  
dataEP.channelName='Random'; %Just to save results
dataEP.numberAntennas=[16 64]; %not used if equalization, but should not be empty 
dataEP.varnoiseCH=0; % If 0, then the channel matrix is perfectly known. In other
%case, this will be the noise variance of the channel matrix
dataEP.numberChannels=100;  %number of random channels. If channel_h is not empty
%this parameter is not used
if isempty(dataEP.channel_h) %Random channels
    dataEP.numberTaps=7; % is the number of taps of the channel to be randomly generated
else
    dataEP.numberTaps=length(dataEP.channel_h);
end

%dataEP.folderSave='ResultsBER';

scenario.flagNorm=1; % If channel is normalized
% Particular Range of SNRdB for algorithms
% Must be within the previous general range

% Methods to simulate
scenario.flagMMSE=1; % block LMMSE [Muranov10]
scenario.flagOptimal=0; % Optimal solution (MAP (BCJR in Equalization))
scenario.flagBEP=0; 
scenario.flagPBEP=1; 
scenario.flagDBEP=1; 
switch dataEP.scenario
    case 1,
        scenario.System='Equalization';
        %Particular approaches to be simulated
        scenario.flagFMMSE=0; % F-LMMSE [Tuchler02]
        scenario.flagPFEP=1; % P-FEP [Santos18]
        scenario.flagDFEP=1; % D-FEP [Santos19]
        scenario.flagSEP=0; % SEP [Santos17]
        scenario.flagPKSEP=0; % P-KSEP [Santos18c]
        scenario.flagDKSEP=0; % D-KSEP [Santos19]
        scenario.flagBPEP=1; % BP-EP [Sun15]
        
        %MIMO Approaches are not simulated
        scenario.flagMMSEAproxGS=0;
        scenario.flagMMSEAproxNeuman=0; % approximated block LMMSE with Neuman
        scenario.flagMMSEAproxPCG=0; % approximated block LMMSE with PCG
        scenario.flagDBEPAproxGS=0; % approximated D-BEP with GS
        scenario.flagDBEPAproxNeuman=0; % approximated D-BEP with Neuman
        scenario.flagDBEPAproxPCG=0; % approximated D-BEP with PCG
        scenario.flagEPICLMMSE=0; % EP IC-LMMSE [Senst11]
                
    case 2
        scenario.System='MIMO';
        scenario.flagMMSEAproxGS=1; % approximated block LMMSE with GS
        scenario.flagMMSEAproxNeuman=0; % approximated block LMMSE with Neuman
        scenario.flagMMSEAproxPCG=0; % approximated block LMMSE with PCG
        scenario.flagDBEPAproxGS=1; % approximated D-BEP with GS
        scenario.flagDBEPAproxNeuman=1; % approximated D-BEP with Neuman
        scenario.flagDBEPAproxPCG=1; % approximated D-BEP with PCG
        scenario.flagEPICLMMSE=1; % EP IC-LMMSE [Senst11]
        
        %Equalization approaches are not simulated
        scenario.flagFMMSE=0; % F-LMMSE [Tuchler02]
        scenario.flagPFEP=0; % P-FEP [Santos18]
        scenario.flagDFEP=0; % D-FEP [Santos19]
        scenario.flagSEP=0; % SEP [Santos17]
        scenario.flagPKSEP=0; % P-KSEP [Santos18c]
        scenario.flagDKSEP=0; % D-KSEP [Santos19]
        scenario.flagBPEP=0; % BP-EP [Sun15]       
        
    otherwise,
        if verbose, disp('Note that no specific approaches have been selected to be simulated'),
%             count=fwrite(fidLog,'No scenario found, EXIT','char');
%             count=fwrite(fidLog,10,'uint8');
%             fclose(fidLog)
        end
%        return
end