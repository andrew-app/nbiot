numTrBlks = 100;        % Number of simulated transport blocks
SNRdB = -35:1:0;      % SNR range in dB
ireps = [0]; 

NPDSCHDataType = 'NotBCCH';

ISF = 0;                % Resource assignment field in DCI (DCI format N1 or N2)
SchedulingInfoSIB1 = 0; % Scheduling information field in MasterInformationBlock-NB (MIB-NB)
IMCS = 4;

enb.NFrame = 0;     % Simulation starting frame number
enb.NSubframe = 0;  % Simulation starting subframe number
enb.NNCellID = 0;   % NB-IoT physical cell ID
enb.NBRefP = 2;     % Number of NRS antenna ports, should be either 1 or 2
enb.OperationMode = 'Inband-DifferentPCI';  % The allowed values are 'Inband-SamePCI', 'Inband-DifferentPCI', 'Guardband' or 'Standalone'
if strcmpi(enb.OperationMode,'Inband-SamePCI')
    enb.CellRefP = enb.NBRefP;     % The allowed values are NBRefP or 4
    enb.NCellID = enb.NNCellID;
elseif strcmpi(enb.OperationMode,'Inband-DifferentPCI')
    enb.CellRefP = 4; % Number of Cell RS antenna ports (Must be equal to NBRefP or 4)      
    enb.NCellID = 1;
end
if (strcmpi(NPDSCHDataType,'BCCHNotSIB1NB') || strcmpi(NPDSCHDataType,'NotBCCH')) && ...
        (strcmpi(enb.OperationMode,'Inband-SamePCI') || strcmpi(enb.OperationMode,'Inband-DifferentPCI'))
    enb.ControlRegionSize = 3;     % The allowed values are 0...13
end

channel = struct;                    % Initialize channel config structure
channel.Seed = 6;                    % Channel seed
channel.NRxAnts = 1;                 % 1 receive antenna
channel.DelayProfile ='EPA';         % Delay profile
channel.DopplerFreq = 5;             % Doppler frequency in Hz
channel.MIMOCorrelation = 'Low';     % Multi-antenna correlation
channel.NTerms = 16;                 % Oscillators used in fading model
channel.ModelType = 'GMEDS';         % Rayleigh fading model type
channel.InitPhase = 'Random';        % Random initial phases
channel.NormalizePathGains = 'On';   % Normalize delay profile power     
channel.NormalizeTxAnts = 'On';      % Normalize for transmit antennas

perfectChannelEstimator = true;

cec.PilotAverage = 'UserDefined';   % Type of pilot symbol averaging
cec.TimeWindow = 1;                 % Time window size in REs
cec.FreqWindow = 25;                % Frequency window size in REs
cec.InterpType = 'Cubic';           % 2D interpolation type
cec.InterpWindow = 'Centered';      % Interpolation window type
cec.InterpWinSize = 3;              % Interpolation window size
cec.Reference = 'NRS';              % Channel estimator reference signal

for repIdx = 1:numel(ireps)
    
    npdschInfo = hNPDSCHInfo;
    npdschInfo.NPDSCHDataType = NPDSCHDataType;
    npdschInfo.ISF = ISF;
    if strcmpi(NPDSCHDataType,'SIB1NB')  % NPDSCH carrying SIB1-NB
        npdschInfo.SchedulingInfoSIB1 = SchedulingInfoSIB1;
    else % NPDSCH not carrying SIB1-NB
        npdschInfo.IRep = ireps(repIdx); % Repetition number field in DCI (DCI format N1 or N2)
        npdschInfo.IMCS = IMCS;          % Modulation and coding scheme field in DCI (DCI format N1 or N2)
        % Verify the inputs of IRep and IMCS
        if isempty(npdschInfo.TBS)
            npdschInfo.TBSTable
            error(['Invalid [ITBS,ISF] (where ITBS=IMCS=' num2str(IMCS)...
                ', ISF=' num2str(ISF)  ') pair, empty TBS is returned, check valid pairs in the above table or 3GPP TS 36.213 table 16.4.1.5.1-1']);
        end
    end
    
    npdsch.NSF = npdschInfo.NSF;
    npdsch.NRep = npdschInfo.NRep;
    npdsch.NPDSCHDataType = NPDSCHDataType;
    npdsch.RNTI = 1;
    
    [~,info] = lteNPDSCHIndices(enb,npdsch);
    rmoutlen = info.G;           % Bit length after rate matching, i.e. codeword length
    trblklen = npdschInfo.TBS;   % Transport block size
    R = (trblklen+24)/rmoutlen;  % DL-SCH channel coding rate, 24 denotes the number of CRC bits
    if R >= 1
        error(['DL-SCH coding rate (' num2str(R) ') larger than or equal to 1 for the configured parameters.']);
    end
    
     %  The NPDSCH repetition pattern for the current configuration is
    %  displayed below
    displayPattern = false;
    % Display NPDSCH repetition pattern
    if displayPattern == true
        npdschInfo.displaySubframePattern;
    end

     % Absolute subframe number at the starting point of the simulation
    NSubframe = enb.NFrame*10+enb.NSubframe;      

    % Initialize BLER and throughput result
    maxThroughput = zeros(length(SNRdB),1);
    simThroughput = zeros(length(SNRdB),1);
    bler = zeros(1,numel(SNRdB));                   

    % The temporary variables 'enb_init' and 'channel_init' are used to create
    % the temporary variable 'enb' and 'channel' within the SNR loop to create
    % independent simulation loops for the 'parfor' loop
    enb_init = enb;
    channel_init = channel;

    %for snrIdx = 1:numel(SNRdB)
    parfor snrIdx = 1:numel(SNRdB)
        
         % Set the random number generator seed depending to the loop variable
        % to ensure independent random streams
        rng(snrIdx,'combRecursive');

        fprintf('\nSimulating %d transport blocks at %gdB SNR\n',numTrBlks,SNRdB(snrIdx));

        enb = enb_init;         % Initialize eNodeB configuration
        channel = channel_init; % Initialize fading channel configuration
        txcw = [];              % Initialize the transmitted codeword
        numBlkErrors = 0;       % Number of transport blocks with errors
        estate = [];            % Initialize NPDSCH encoder state
        dstate = [];            % Initialize NPDSCH decoder state
        lastOffset = 0;         % Initialize overall frame timing offset
        offset = 0;             % Initialize frame timing offset
        subframeGrid = lteNBResourceGrid(enb); % Initialize the subframe grid

        subframeIdx = NSubframe;
        numRxTrBlks = 0;
        while (numRxTrBlks < numTrBlks)

            % Set current subframe and frame numbers  
            enb.NSubframe = mod(subframeIdx,10);
            enb.NFrame = floor((subframeIdx)/10);
            
            % Generate the NPSS symbols and indices
            npssSymbols = lteNPSS(enb);
            npssIndices = lteNPSSIndices(enb);
            % Map the symbols to the subframe grid
            subframeGrid(npssIndices) = npssSymbols;
            
            % Generate the NSSS symbols and indices
            nsssSymbols = lteNSSS(enb);
            nsssIndices = lteNSSSIndices(enb);
            % Map the symbols to the subframe grid
            subframeGrid(nsssIndices) = nsssSymbols;
            
            % Establish if either NPSS or NSSS is transmitted and if so,
            % do not transmit NPDSCH in this subframe
            isDataSubframe = isempty(npssSymbols) && isempty(nsssSymbols);

            % Create a new transport block and encode it when the
            % transmitted codeword is empty. The receiver sets the codeword
            % to empty to signal that all subframes in a bundle have been
            % received (it is also empty before the first transmission)
            if isempty(txcw)
                txTrBlk = randi([0 1],trblklen,1);
                txcw = lteNDLSCH(rmoutlen,txTrBlk);
            end

            if (isDataSubframe)
                % Generate NPDSCH symbols and indices for a subframe
                [txNpdschSymbols,estate] = lteNPDSCH(enb,npdsch,txcw,estate);
                npdschIndices = lteNPDSCHIndices(enb,npdsch);
                % Map the symbols to the subframe grid
                subframeGrid(npdschIndices) = txNpdschSymbols;
                % Generate the NRS symbols and indices
                nrsSymbols = lteNRS(enb);
                nrsIndices = lteNRSIndices(enb);
                % Map the symbols to the subframe grid 
                subframeGrid(nrsIndices) = nrsSymbols;
            end

            % Perform OFDM modulation to generate the time domain waveform
            [txWaveform,ofdmInfo] = nbOFDMModulate(enb,subframeGrid);

            % Add 25 sample padding. This is to cover the range of delays
            % expected from channel modeling (a combination of
            % implementation delay and channel delay spread)
            txWaveform =  [txWaveform; zeros(25, enb.NBRefP)]; %#ok<AGROW>

            % Initialize channel time for each subframe
            channel.InitTime = subframeIdx/1000;

            % Pass data through channel model
            channel.SamplingRate = ofdmInfo.SamplingRate;
            [rxWaveform,fadingInfo] = lteFadingChannel(channel, txWaveform);

            % Calculate noise gain including compensation for downlink power
            % allocation
            SNR = 10^(SNRdB(snrIdx)/20);

            % Normalize noise power to take account of sampling rate, which
            % is a function of the IFFT size used in OFDM modulation, and
            % the number of antennas
            N0 = 1/(sqrt(2.0*enb.NBRefP*double(ofdmInfo.Nfft))*SNR);

            % Create additive white Gaussian noise
            noise = N0*complex(randn(size(rxWaveform)), ...
                                randn(size(rxWaveform)));

            % Add AWGN to the received time domain waveform        
            rxWaveform = rxWaveform + noise;
            
            if(perfectChannelEstimator)
                offset = hPerfectTimingEstimate(fadingInfo);
            else
                % In this example, the subframe offset calculation relies
                % on NPSS present in subframe 5, so we need to pad the
                % subframes before it so that the frame offset returned by
                % lteNBDLFrameOffset is the offset for subframe 5
                sfTsamples = ofdmInfo.SamplingRate*1e-3;
                if (enb.NSubframe==5) 
                    padding = zeros([sfTsamples*5,size(rxWaveform,2)]);
                    offset = lteNBDLFrameOffset(enb, [padding; rxWaveform]);
                    if (offset > 25) || (offset < 0)
                        offset = lastOffset;
                    end
                    lastOffset = offset;
                end
            end

            % Synchronize the received waveform
            rxWaveform = rxWaveform(1+offset:end, :);

            % Perform OFDM demodulation on the received data to recreate the
            % resource grid
            rxSubframe = nbOFDMDemodulate(enb,rxWaveform);
            
            % Channel estimation
            if(perfectChannelEstimator) 
                % Perfect channel estimation
                estChannelGrid = nbDLPerfectChannelEstimate(enb, channel, offset);
                noiseGrid = nbOFDMDemodulate(enb, noise(1+offset:end ,:));
                noiseEst = var(noiseGrid(:));
            else

                [estChannelGrid, noiseEst] = lteDLChannelEstimate( ...
                enb, cec, rxSubframe);
            end

            if (isDataSubframe)
                % Get NPDSCH indices
                npdschIndices = lteNPDSCHIndices(enb, npdsch);

                % Get PDSCH resource elements from the received subframe. Scale the
                % received subframe by the PDSCH power factor Rho. The PDSCH is
                % scaled by this amount, while the cell reference symbols used for
                % channel estimation (used in the PDSCH decoding stage) are not.
                [rxNpdschSymbols, npdschHest] = lteExtractResources(npdschIndices, ...
                    rxSubframe, estChannelGrid);

                % Decode NPDSCH
                [rxcw,dstate,symbols] = lteNPDSCHDecode(...
                                     enb, npdsch, rxNpdschSymbols, npdschHest, noiseEst,dstate);

                % Decode the transport block when all the subframes in a bundle
                % have been received
                if dstate.EndOfTx
                   [trblkout,blkerr] = lteNDLSCHDecode(trblklen,rxcw);
                   numBlkErrors = numBlkErrors + blkerr;
                   numRxTrBlks = numRxTrBlks + 1;
                   % Re-initialize to enable the transmission of a new transport block
                   txcw = [];
                end
            end

            subframeIdx = subframeIdx + 1;
            
        end

        % Calculate the block error rate
        bler(snrIdx) = numBlkErrors/numTrBlks;
        fprintf('NPDSCH BLER = %.4f \n',bler(snrIdx));
        % Calculate the maximum and simulated throughput
        maxThroughput(snrIdx) = trblklen*numTrBlks; % Max possible throughput
        simThroughput(snrIdx) = trblklen*(numTrBlks-numBlkErrors);  % Simulated throughput
        fprintf('NPDSCH Throughput(%%) = %.4f %%\n',simThroughput(snrIdx)*100/maxThroughput(snrIdx));

    end
    
    if repIdx == 1
        fh = figure;
        grid on;
        hold on;
        xlabel('SNR (dB)');
        ylabel('BLER');
        legendstr = {['NRep = ' num2str(npdsch.NRep)]};
    else
        legendstr = [legendstr ['NRep = ' num2str(npdsch.NRep)]]; %#ok<AGROW>
    end
    figure(fh);
    plot(SNRdB, bler, '-o');


end
% Set figure title
if strcmpi(NPDSCHDataType,'SIB1NB')
    npdsch.NSF = 8;
end
title([' ' char(npdsch.NPDSCHDataType) ': TBS=' num2str(trblklen)...
    '; NSF=' num2str(npdsch.NSF) '; ' num2str(enb_init.NBRefP) ' NRS port(s)' ]);
legend(legendstr);

%% Local functions 

% NB-IoT DL OFDM Modulator
function [waveform,info] = nbOFDMModulate(enb,grid)
    % Apply default window size according to TS 36.104 Table E.5.1-1a
    if(~isfield(enb,'Windowing'))
        enb.Windowing = 6;
    end
    % Use NB-IoT SC-FDMA to get the 1/2 subcarrier shift on the OFDM modulation
    enb.NBULSubcarrierSpacing = '15kHz'; 
    [waveform,info] = lteSCFDMAModulate(enb,grid);
end

% NB-IoT DL OFDM Demodulator
function grid = nbOFDMDemodulate(enb,rxWaveform)
    % Use NB-IoT SC-FDMA to get the 1/2 subcarrier shift on the OFDM modulation
    enb.NBULSubcarrierSpacing = '15kHz'; 
    grid = lteSCFDMADemodulate(enb,rxWaveform,0.55); % CP fraction of 0.55
end

% NB-IoT DL Perfect Channel Estimator
function H = nbDLPerfectChannelEstimate(enb,channel,timefreqoffset)
    % Reconfigure NB-IoT UL perfect channel estimator to perform DL perfect
    % channel estimation
    enb.NBULSubcarrierSpacing = '15kHz'; 
    enb.NTxAnts = enb.NBRefP;
    enb.TotSlots = 2; 
    H = lteULPerfectChannelEstimate(enb, channel,timefreqoffset);
end
