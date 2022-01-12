clear variables
clc;
single_test_run = 0; % single SNR (for debugg or dNN dataset)
print_logs = 1; % prints some info regarding the simulation status
save_workspace = 1; % after finishing the simulation, save the current workspace into a temporary file. (The data would be overwritten)

% General parameters
snr = [-5:5:30]
its = 10; % number of simulation iterations for each snr (must be >= 2)
genie_sync = 1; % synchronization flag 1 (true), 0 (false)
nT = 1; % number of Tx antennas
nR = 1; % number of Rx antennas
mu = 4; % Modulation order 2^mu-QAM
noblk = 120; % number of blocks (OFDM symbols)
midamble_spacing = 4; % time domain midamble symbol spacing for 11.bd 
no_midambles = floor(noblk/(midamble_spacing-1)); % number of midambles

% OFDM parameters
N = 64; % FFT size
N_data = 48; % Number of data subcarriers
Np = 4; % Number of pilot subcarriers
Fs = 10e6; % Sampling Frequency
ShortCP = 1.6e-6; % Short CP in seconds
LongCP = 3.2e-6; % Long CP in seconds
preamble_T1 = 6.4e-6; % Duration of T1 of preamble in seconds
Ncp = ShortCP*Fs; % number of CP samples for short CP
Ncp_long = LongCP*Fs; % number of CP samples for long CP
N_Ncp = N + Ncp; % symbol length with CP
T1_length = preamble_T1*Fs; % number of samples in T1
data_subcar_indexes_0 = [-26:-22 -20:-8 -6:-1 1:6 8:20 22:26].'; % set of active data subcarriers
data_subcar_indexes = 7+26+data_subcar_indexes_0; 
pilot_subcar_indexes_0 = [-21 -7 7 21]; % set of active pilot subcarriers
pilot_subcar_indexes = 7+26+pilot_subcar_indexes_0;


% Unique Word parameters
PrimeNr = 5;  % Prime number for Polyphase sequence
UW_length = Ncp*(nT+1); % Length in samples
N_Nuw = UW_length + N; % one block/symbol duration of payload-UW
frame_length_UW = noblk*N_Nuw;
UW_length_withcp = Ncp + PrimeNr^2;
N_Nuwcp = UW_length_withcp + N;

% Channel parameters
path_delays = [1:8];  % Path Delays (At the moment integer multiple of 1/Fs)
numPaths = length(path_delays); % Number of Channel Paths
AveGaindB = linspace(0,-20,numPaths); % Gain of the channel paths in dB
numTaps = numPaths; % Number of non-zero channel taps 
n_H = numPaths + round((UW_length-numPaths)/2); % (used for Wiener filtering) sample index where CE is at its best
n_Hcpuw = Ncp + round(PrimeNr^2/2);
n_H_11bd = Ncp + round(N/2); % (used for Wiener filtering) sample index where CE is at its best
P = get_pdp(nT,nR,numPaths,AveGaindB,'exp decaying'); % Power Delay Profile
preamble.R_hh = kron(eye(nT),diag(P(:,1,1))); % Covariance mtx for channel estimation
chan_typ = 'time variant'; 
fd = 1600; % max Doppler Shift


if single_test_run
   snr = [35]
  
   fprintf('TRY RUN with snr %d dB, fd = %4.2f \n',snr(1),fd)
   fprintf('TRY RUN with snr %d dB, fd = %4.2f \n',snr(1),fd)
   fprintf('TRY RUN with snr %d dB, fd = %4.2f \n',snr(1),fd)
end


% Coding parameters
trellis= poly2trellis(7,[133 171]); % Construct a trellis for Turbo encoder
n = log2(trellis.numOutputSymbols);
numTails = log2(trellis.numStates)*n;
N_outBits = N_data*mu; % Number of output codeword bits
interleaver_size = (N_outBits - 2*numTails)/(2*n-1); % Size of code interleaver
coderate = interleaver_size/N_outBits; % Code rate
clear n numTails

% Add the current variables into a structure
globpara=[who,who].'; 
config=struct(globpara{:}); 
eval(structvars(config,0).')  % All the current variables are in config struct

% Get the UW sequence
UW_inTime_1 = zeros(UW_length,nT);
UW_inTime_withcp = zeros(PrimeNr^2,nT);
for iT=1:nT
%    UW_inTime_1(:,iT) = get_uw_sequence( UW_length, iT, 'Polyphase', PrimeNr); 
   UW_inTime_1(:,iT) = get_uw_sequence( UW_length, iT, 'Zadoff-Chu'); 
   
   UW_inTime_withcp(:,iT) = get_uw_sequence( PrimeNr^2, iT, 'Polyphase', PrimeNr); 
end
config.UW_inFreq = fft_u(UW_inTime_withcp);
UW_inTime_withcp = do_insert_CP(UW_inTime_withcp,Ncp);
config.UW_inTime_withcp = UW_inTime_withcp;

% Get Preamble and Pilots for CP-OFDM Channel Estimation
preambleT1_inTime = zeros(T1_length+1,nT);
pilots_inTime = zeros(N,nT);
midamble_inTime = zeros(T1_length+1,nT,no_midambles);
for iT=1:nT
   % For both 802.11p and 802.11bd
   preambleT1_inTime(:,iT) = zadoffChuSeq(iT,T1_length+1);
   pilots_inTime = get_80211p_pilots(config);
   % For 802.11bd
   for midamble_id = 1:no_midambles
       midamble_inTime(:,iT,midamble_id) = zadoffChuSeq(1,T1_length+1); 
   end
   
  
end
preambleT1_inTime(end,:,:)=[];
preambleT1_inFreq = fft_u(preambleT1_inTime);
preamble_inTime_withCP = do_insert_CP(cat(1,preambleT1_inTime,preambleT1_inTime),Ncp_long);
pilots_inTime_withCP = do_insert_CP(pilots_inTime,Ncp);
midamble_inTime(end,:,:)=[];
config.midamble_inFreq = fft_u(midamble_inTime);
midamble_withCP = do_insert_CP(midamble_inTime,Ncp);
 
% Construct Toeplitz observation matrix for channel estimation
config = get_obs_mtx(config,UW_inTime_1,UW_length,numTaps,nT);

ber_CP = zeros(length(snr),its);
ber_11bd_noW = zeros(length(snr),its);
ber_11bd = zeros(length(snr),its);
ber_UW = zeros(length(snr),its);
ber_UW_noW = zeros(length(snr),its);
ber_UWcp = zeros(length(snr),its);
ber_UWcp_noW = zeros(length(snr),its);
mse_CP = zeros(length(snr),its);
mse_UW = zeros(length(snr),its);
mse_UW_noW = zeros(length(snr),its);
mse_UWcp = zeros(length(snr),its);
mse_UWcp_noW = zeros(length(snr),its);
mse_11bd_noW = zeros(length(snr),its);
mse_11bd = zeros(length(snr),its);

Raw_DNN_dataset.its = its;
Raw_DNN_dataset.N = N;
Raw_DNN_dataset.noblk = noblk;
Raw_DNN_dataset.snr = snr;


start_time = datestr(now,'dd.mm.yyyy at HH:MM:SS');
for si=1:length(snr)
    
    Raw_DNN_dataset.h_genie_UW_dnn = [];
    Raw_DNN_dataset.h_genie_CP_11bd_dnn = [];
    Raw_DNN_dataset.h_genie_UWcp_dnn = [];
    Raw_DNN_dataset.H_hat_out_noW_dnn = [];
    Raw_DNN_dataset.H_hat_out_dnn = [];
    Raw_DNN_dataset.H_hat_out_11bd_noW_dnn = [];
    Raw_DNN_dataset.H_hat_out_11bd_dnn = [];
    Raw_DNN_dataset.H_hat_out_uwcp_noW_dnn = [];
    Raw_DNN_dataset.H_hat_out_uwcp_dnn = [];
    
    for i = 1:its
        
        if print_logs
            if mod(i,round(its/4)) == 0
                fprintf('at %s, running at iteration %d.\n',datestr(now,'HH:MM:SS'),i)
            elseif i == 1
                fprintf('Simulating SNR %d dB at max Doppler shift %d Hz. \nTime: %s\n',snr(si),fd,datestr(now,'HH:MM:SS'))
            end
        end
    
        % Generate the channel(s) filter
        h = get_mimo_channel(P(1:numPaths,:,:),Fs,numPaths,path_delays,nT,nR,chan_typ,fd);
        
        % Create some dummy data for the start and end of the block
        dummy_length = 200;
        dummy_end = (randn(dummy_length, nT)+1i*randn(dummy_length, nT))/sqrt(2);
        dummy_start = (randn(dummy_length, nT)+1i*randn(dummy_length, nT))/sqrt(2);
        
%% CP-based frames 

        % Get the OFDM payloads without CP
        [payload_inTime, config] = get_payload_ofdm(config, noblk, 'coded');
        
        % CP Insertion for payload with CP
        payload_withCP = do_insert_CP(payload_inTime,Ncp);
        
        % Pilot Insertion for payload with CP
        payload_withCP = payload_withCP + pilots_inTime_withCP;   
        
        % CP-based Tx frame: concatenate all CP-Payload blocks (802.11p)
        tx_signal_withCP = [dummy_start; preamble_inTime_withCP];
        for blk = 1:noblk
            tx_signal_withCP = [tx_signal_withCP; payload_withCP(:,:,blk)];
        end
        tx_signal_withCP = [tx_signal_withCP; dummy_end];       
                
        % filter the signals with the generated channels
        data_middle_points_11p = (dummy_length+1) + length(preamble_inTime_withCP) + Ncp + N_Ncp*(0:noblk-1) + N/2; % in order to extract the genie-aided CSI
        [xch_withCP, h_genie_CP] = do_mimo_channel('Filter time variant',tx_signal_withCP, h, numPaths, nT, nR, length(tx_signal_withCP), noblk, data_middle_points_11p);
        
        % AWGN process to the channel filtered signal
        w_cp = awgn(complex(zeros(length(xch_withCP(:,1)),nT)),snr(si)); % noise is being separately generated just to keep track of it while debugging
        rx_signal_withCP = xch_withCP + w_cp; 
               
        % Synchronization
        if genie_sync
            preamble_starts = dummy_length+1;
            CP_payload_starts = (preamble_starts + Ncp_long + 2*T1_length + Ncp):N_Ncp:(preamble_starts + Ncp_long + 2*T1_length + Ncp)+(noblk-1)*N_Ncp;
        else
            % Synchronization errors can be added here
        end
        
        % Primary channel estimation based on preamble (CP-based packets)
        rx_preamble_CP = rx_signal_withCP(preamble_starts+Ncp_long:preamble_starts+Ncp_long+2*T1_length-1,:);
        H_hat_in_cpofdm = do_fft_channel_estimation(config,rx_preamble_CP,preambleT1_inFreq,snr(si),'true','low complexity approx. LMMSE');
        
        % Exctract the N-point data + cp
        y_cpofdm = zeros(N,nR,noblk);        
        for blk = 1:noblk
            y_cpofdm(:,:,blk) = rx_signal_withCP(CP_payload_starts(blk):CP_payload_starts(blk)+N-1,:);
            %     x2(:,:,blk) = tx_signal_part2(data_starts(blk):data_starts(blk)+N_Nuw-1,:); % For debugging
        end                
        
        % CP-OFDM equalization based on preamble
        d_hat_cpofdm = cell(noblk,1);
        H_hat_out_cpofdm = cell(noblk,1);
        Rx_matrices_cpofdm.R_dd_hat = cell(noblk,1);
        for blk = 1:noblk
            Y_cpofdm = fft_u(y_cpofdm(:,:,blk));
            Y_cpofdm_v = Y_cpofdm(:);
            [ d_hat_cpofdm{blk} , ~, ~, R_dd_hat_cpofdm, ~] = do_mimo_mmse_equalization_cpofdm( config, snr(si), H_hat_in_cpofdm, zeros(numPaths), Y_cpofdm_v);
            Rx_matrices_cpofdm.R_dd_hat{blk} = diag(R_dd_hat_cpofdm);
            H_hat_out_cpofdm{blk} = H_hat_in_cpofdm;
        end
                
        % Soft-out demodulation and decoding
        [receivedBits_CP, ~, ~] = do_soft_decode(config, d_hat_cpofdm, Rx_matrices_cpofdm);

%% IEEE 802.11bd frames

        % CP-based Tx frame: concatenate all CP-Payload blocks (802.11bd)
        tx_signal_withCP_11bd = [dummy_start; preamble_inTime_withCP];
        for midamble_id = 0:(no_midambles-1)
            tx_signal_withCP_11bd = [tx_signal_withCP_11bd; payload_withCP(:,:,midamble_id*(midamble_spacing-1) + 1); ...
               payload_withCP(:,:,midamble_id*(midamble_spacing-1)+2); payload_withCP(:,:,midamble_id*(midamble_spacing-1)+3); ...
               midamble_withCP(:,:,midamble_id+1)];
        end
        if (mod(noblk,3) == 1)
            tx_signal_withCP_11bd = [tx_signal_withCP_11bd; payload_withCP(:,:,end)];
        elseif (mod(noblk,3) == 2)
            tx_signal_withCP_11bd = [tx_signal_withCP_11bd; payload_withCP(:,:,end-1); payload_withCP(:,:,end)];
        end
        tx_signal_withCP_11bd = [tx_signal_withCP_11bd; dummy_end];
        
        % filter the signals with the generated channels
        data_middle_points_11bd = (dummy_length+1) + length(preamble_inTime_withCP) + Ncp + N_Ncp*(0:noblk+no_midambles-1) + N/2; % in order to extract the genie-aided CSI
        data_middle_points_11bd(4:4:end) = []; % remove the midamble indexes
        [xch_withCP_11bd, h_genie_CP_11bd] = do_mimo_channel('Filter time variant',tx_signal_withCP_11bd, h, numPaths, nT, nR, length(tx_signal_withCP_11bd), noblk, data_middle_points_11bd);

        % AWGN process to the channel filtered signal
        w_cp_11bd = awgn(complex(zeros(length(xch_withCP_11bd(:,1)),nT)),snr(si));
        rx_signal_withCP_11bd = xch_withCP_11bd + w_cp_11bd;
               
        % Synchronization
        if genie_sync
            preamble_starts = dummy_length+1;
            CP_symbol_starts_11bd = (preamble_starts + Ncp_long + 2*T1_length + Ncp):N_Ncp:(preamble_starts + Ncp_long + 2*T1_length + Ncp)+(noblk+no_midambles-1)*N_Ncp;          
        else
            % Synchronization errors can be added here
        end
        
        % Primary channel estimation based on preamble (CP-based packets)
        rx_preamble_CP_11bd = rx_signal_withCP_11bd(preamble_starts+Ncp_long:preamble_starts+Ncp_long+2*T1_length-1,:);
        [h_hat_primary_11bd, R_hh_err_11bd] = do_fft_channel_estimation(config,rx_preamble_CP_11bd,preambleT1_inFreq,snr(si),'true','exact LMMSE');
%         h_hat_primary_11bd = ifft(H_hat_in_cpofdm_11bd); h_hat_primary_11bd(numPaths+1:end) = [];
        
        % Exctract the N-point data + cp
        y_cpofdm_11bd = zeros(N,nR,noblk+no_midambles);        
        for blk = 1:(noblk+no_midambles)
            y_cpofdm_11bd(:,:,blk) = rx_signal_withCP_11bd(CP_symbol_starts_11bd(blk):CP_symbol_starts_11bd(blk)+N-1,:);
            %     x2(:,:,blk) = tx_signal_part2(data_starts(blk):data_starts(blk)+N_Nuw-1,:); % For debugging
        end
        
        % Joint Channel Estimation and Equalization
        wiener_flag = 0;
        [ d_hat_11bd_noW, H_hat_out_11bd_noW, Rx_matrices_11bd_noW ] = do_joint_ce_eq_11bd(config, snr(si), y_cpofdm_11bd, h_hat_primary_11bd, R_hh_err_11bd, noblk, wiener_flag);
        wiener_flag = 2; % with prediction/extrapolation after the last midamble
        [ d_hat_11bd, H_hat_out_11bd, Rx_matrices_11bd ] = do_joint_ce_eq_11bd(config, snr(si), y_cpofdm_11bd, h_hat_primary_11bd, R_hh_err_11bd, noblk, wiener_flag);
        
        % Soft-out demodulation and decoding
        [receivedBits_11bd_noW, ~, ~] = do_soft_decode(config, d_hat_11bd_noW, Rx_matrices_11bd_noW);
        [receivedBits_11bd, ~, ~] = do_soft_decode(config, d_hat_11bd, Rx_matrices_11bd);
        
%% Unique-Word based frames

        % UW Insertion for payload with UW
        payload_withUW = do_insert_UW(payload_inTime,UW_inTime_1);
        
        % UW-based Tx frame: concatenate all UW-Payoad blocks
        tx_signal_withUW = [dummy_start; preamble_inTime_withCP];
        for blk = 1:noblk
            tx_signal_withUW = [tx_signal_withUW; payload_withUW(:,:,blk)];
        end
        tx_signal_withUW = [tx_signal_withUW; UW_inTime_1; dummy_end];
        
        % filter the signals with the generated channels
        data_middle_points_uw = (dummy_length+1) + length(preamble_inTime_withCP) + N_Nuw*(0:noblk-1) + N/2; % in order to extract the genie-aided CSI
        [xch_withUW, h_genie_UW] = do_mimo_channel('Filter time variant',tx_signal_withUW, h, numPaths, nT, nR, length(tx_signal_withUW), noblk, data_middle_points_uw);
        
        % AWGN process to the channel filtered signal
        w_uw = awgn(complex(zeros(length(xch_withUW(:,1)),nT)),snr(si));
        rx_signal_withUW = xch_withUW + w_uw;
               
        % Synchronization
        if genie_sync
            preamble_starts = dummy_length+1;
            %UW_payload_starts = (preamble_starts + 2*length(UW_inTime_1)):N_Nuw:(preamble_starts + 2*length(UW_inTime_1))+(noblk-1)*N_Nuw;            
            UW_payload_starts = (preamble_starts + Ncp_long + 2*T1_length):N_Nuw:(preamble_starts + Ncp_long + 2*T1_length)+(noblk-1)*N_Nuw;            
            
        else
            % Synchronization errors can be added here
        end
        
        % Exctract the N-point data + uw
        y = zeros(N_Nuw,nR,noblk);
        % x2 = zeros(N,nR,noblk);
        for blk = 1:noblk
            y(:,:,blk) = rx_signal_withUW(UW_payload_starts(blk):UW_payload_starts(blk)+N_Nuw-1,:);
            %     x2(:,:,blk) = tx_signal_part2(data_starts(blk):data_starts(blk)+N_Nuw-1,:); % For debugging
        end
        
%         % Primary channel estimation based on preamble (UW-based packets)
%         rx_preamble_UW = rx_signal_withUW(preamble_starts:preamble_starts+2*length(UW_inTime_1)-1,:);
%         h_hat_primary = zeros(numTaps,nT,nR); % Primary channel estimation
%         R_hh_err = zeros(nT*numTaps,nT*numTaps,nR); % channel estimation covariance matrix
%         for iR = 1:nR
%             % Primary channel estimation
%             [h_hat_primary(:,:,iR), R_hh_err(:,:,iR)] = do_mimo_uw_channel_estimation( config, snr(si), rx_preamble_UW(:,iR), 'primary_toeplitz');
%         end
        
        % Primary channel estimation based on preamble (CP-based packets)
        rx_preamble_UW = rx_signal_withUW(preamble_starts+Ncp_long:preamble_starts+Ncp_long+2*T1_length-1,:);
        [h_hat_primary_UW, R_hh_err_UW] = do_fft_channel_estimation(config,rx_preamble_UW,preambleT1_inFreq,snr(si),'true','exact LMMSE');

        
        % Joint Channel Estimation and Equalization
        wiener_flag = 0;
        [ d_hat_noW, H_hat_out_noW, Rx_matrices_noW ] = do_joint_ce_eq_uwofdm(config, snr(si), y, h_hat_primary_UW, R_hh_err_UW, noblk, wiener_flag);
        wiener_flag = 1;
        [ d_hat, H_hat_out, Rx_matrices ] = do_joint_ce_eq_uwofdm(config, snr(si), y, h_hat_primary_UW, R_hh_err_UW, noblk, wiener_flag);
                              
               
        % Soft-out demodulation and decoding
        [receivedBits_UW_noW, ~, ~] = do_soft_decode(config, d_hat_noW, Rx_matrices_noW); % without Wiener filter
        [receivedBits_UW, ~, ~] = do_soft_decode(config, d_hat, Rx_matrices); % with Wiener filter


%% UW based frame with CP-UW
        
        % UW Insertion for payload with UW (UW includes a CP)
        payload_withUWcp = do_insert_UW(payload_inTime,UW_inTime_withcp);
        
        % UW-based Tx frame: concatenate all UW-Payoad blocks
        tx_signal_withUWcp = [dummy_start; preamble_inTime_withCP];
        for blk = 1:noblk
            tx_signal_withUWcp = [tx_signal_withUWcp; payload_withUWcp(:,:,blk)];
        end
        tx_signal_withUWcp = [tx_signal_withUWcp; UW_inTime_withcp; dummy_end];
        
        % filter the signals with the generated channels
        data_middle_points_uwcp = (dummy_length+1) + length(preamble_inTime_withCP) + N_Nuwcp*(0:noblk-1) + N/2; % in order to extract the genie-aided CSI
        [xch_withUWcp, h_genie_UWcp] = do_mimo_channel('Filter time variant',tx_signal_withUWcp, h, numPaths, nT, nR, length(tx_signal_withUWcp), noblk, data_middle_points_uwcp);
        
        % AWGN process to the channel filtered signal
        w_uwcp = awgn(complex(zeros(length(xch_withUWcp(:,1)),nT)),snr(si));
        rx_signal_withUWcp = xch_withUWcp + w_uwcp;
               
        % Synchronization
        if genie_sync
            preamble_starts = dummy_length+1;
            %UW_payload_starts = (preamble_starts + 2*length(UW_inTime_1)):N_Nuw:(preamble_starts + 2*length(UW_inTime_1))+(noblk-1)*N_Nuw;            
            UWcp_payload_starts = (preamble_starts + Ncp_long + 2*T1_length):N_Nuwcp:(preamble_starts + Ncp_long + 2*T1_length)+(noblk-1)*N_Nuwcp;            
            
        else
            % Synchronization errors can be added here
        end
        
        % Exctract the N-point data + uw
        y_uwcp = zeros(N_Nuwcp,nR,noblk);
        % x2 = zeros(N,nR,noblk);
        for blk = 1:noblk
            y_uwcp(:,:,blk) = rx_signal_withUWcp(UWcp_payload_starts(blk):UWcp_payload_starts(blk)+N_Nuwcp-1,:);
            %     x2(:,:,blk) = tx_signal_part2(data_starts(blk):data_starts(blk)+N_Nuw-1,:); % For debugging
        end
        
        % Primary channel estimation based on preamble (CP-based packets)
        rx_preamble_UWcp = rx_signal_withUWcp(preamble_starts+Ncp_long:preamble_starts+Ncp_long+2*T1_length-1,:);
        [h_hat_primary_UWcp, R_hh_err_UWcp] = do_fft_channel_estimation(config,rx_preamble_UWcp,preambleT1_inFreq,snr(si),'true','exact LMMSE');

        
        % Joint Channel Estimation and Equalization
        wiener_flag = 0;
        [ d_hat_uwcp_noW, H_hat_out_uwcp_noW, Rx_matrices_uwcp_noW ] = do_joint_ce_eq_cpuwofdm(config, snr(si), y_uwcp, h_hat_primary_UWcp, R_hh_err_UWcp, noblk, wiener_flag);
        wiener_flag = 1;
        [ d_hat_uwcp, H_hat_out_uwcp, Rx_matrices_uwcp ] = do_joint_ce_eq_cpuwofdm(config, snr(si), y_uwcp, h_hat_primary_UWcp, R_hh_err_UWcp, noblk, wiener_flag);

        
        % Soft-out demodulation and decoding
        [receivedBits_UWcp_noW, ~, ~] = do_soft_decode(config, d_hat_uwcp_noW, Rx_matrices_uwcp_noW); % without Wiener filter
        [receivedBits_UWcp, ~, ~] = do_soft_decode(config, d_hat_uwcp, Rx_matrices_uwcp); % with Wiener filter

        
        
%% Performance evaluations

        % Calculate the channel estimation MSE     
        mse_UW_noW(si,i) = calc_ce_mse(config,h_genie_UW,H_hat_out_noW);
        mse_UW(si,i) = calc_ce_mse(config,h_genie_UW,H_hat_out);
        mse_CP(si,i) = calc_ce_mse(config,h_genie_CP,H_hat_out_cpofdm);
        mse_11bd_noW(si,i) = calc_ce_mse(config,h_genie_CP_11bd,H_hat_out_11bd_noW);
        mse_11bd(si,i) = calc_ce_mse(config,h_genie_CP_11bd,H_hat_out_11bd);
        mse_UWcp_noW(si,i) = calc_ce_mse(config,h_genie_UWcp,H_hat_out_uwcp_noW);
        mse_UWcp(si,i) = calc_ce_mse(config,h_genie_UWcp,H_hat_out_uwcp);
        
        ber_UW_noW(si,i) = calc_ber(config,receivedBits_UW_noW);
        ber_UW(si,i) = calc_ber(config,receivedBits_UW);
        ber_UWcp_noW(si,i) = calc_ber(config,receivedBits_UWcp_noW);
        ber_UWcp(si,i) = calc_ber(config,receivedBits_UWcp);
        ber_CP(si,i) = calc_ber(config,receivedBits_CP);
        ber_11bd_noW(si,i) = calc_ber(config,receivedBits_11bd_noW);
        ber_11bd(si,i) = calc_ber(config,receivedBits_11bd);
        
        
        % Store the genie channels for DNN (target vectors)
        Raw_DNN_dataset.h_genie_UW_dnn = cat(2,Raw_DNN_dataset.h_genie_UW_dnn,h_genie_UW); clear h_genie_UW;
        Raw_DNN_dataset.h_genie_CP_11bd_dnn = cat(2,Raw_DNN_dataset.h_genie_CP_11bd_dnn,h_genie_CP_11bd); clear h_genie_CP_11bd;
        Raw_DNN_dataset.h_genie_UWcp_dnn = cat(2,Raw_DNN_dataset.h_genie_UWcp_dnn,h_genie_UWcp); clear h_genie_UWcp;
        
        % Store the estimation results for DNN (observations)
        Raw_DNN_dataset.H_hat_out_noW_dnn = cat(2, Raw_DNN_dataset.H_hat_out_noW_dnn, H_hat_out_noW); clear H_hat_out_noW;
        Raw_DNN_dataset.H_hat_out_dnn = cat(2, Raw_DNN_dataset.H_hat_out_dnn, H_hat_out); clear H_hat_out;
        Raw_DNN_dataset.H_hat_out_11bd_noW_dnn = cat(2, Raw_DNN_dataset.H_hat_out_11bd_noW_dnn, H_hat_out_11bd_noW); clear H_hat_out_11bd_noW;
        Raw_DNN_dataset.H_hat_out_11bd_dnn = cat(2, Raw_DNN_dataset.H_hat_out_11bd_dnn, H_hat_out_11bd); clear H_hat_out_11bd;
        Raw_DNN_dataset.H_hat_out_uwcp_noW_dnn = cat(2, Raw_DNN_dataset.H_hat_out_uwcp_noW_dnn, H_hat_out_uwcp_noW); clear H_hat_out_uwcp_noW;
        Raw_DNN_dataset.H_hat_out_uwcp_dnn = cat(2, Raw_DNN_dataset.H_hat_out_uwcp_dnn, H_hat_out_uwcp); clear H_hat_out_uwcp;
        Raw_DNN_dataset.y_uwcp{i} = y_uwcp;         
        Raw_DNN_dataset.y_cpofdm_11bd{i} = y_cpofdm_11bd;             
        Raw_DNN_dataset.Rx_matrices_uwcp_noW{i} = Rx_matrices_uwcp_noW;
        Raw_DNN_dataset.Rx_matrices_11bd_noW{i} = Rx_matrices_11bd_noW;
        

    end
    if print_logs
        fprintf('Simulation of SNR %d ended.\n',snr(si))
        if si==length(snr)
            fprintf('Time: %s\n',datestr(now,'HH:MM:SS'))
        end
    end
    if save_workspace
        save(['Raw_DNN_dataset_fd_' num2str(fd) '_snr_', num2str(si)], '-struct','Raw_DNN_dataset')
        fprintf('\nSimulation data stored in \n%s\\Raw_DNN_dataset_fd_%d_snr_%d.mat\n',pwd,fd,si)
    end
end
end_time = datestr(now,'dd.mm.yyyy at HH:MM:SS');

% SNR per data symbol (Taking into account the overheads)
snr_UW = snr +(UW_length + N)/N;
snr_UWcp = snr +(UW_length_withcp + N)/N;
snr_11p = snr + (Ncp + N)/N;
snr_11bd = snr + ((noblk+no_midambles)*(Ncp+N))/(noblk*N); 

% color='rmbc';
markers = '+o*xv^p';

figure;
semilogy(snr_UW,abs(mean(mse_UW,2)),'Marker',markers(1))
hold on 
semilogy(snr_UW,abs(mean(mse_UW_noW,2)),'Marker',markers(2))
semilogy(snr_11p,abs(mean(mse_CP,2)),'Marker',markers(3))
semilogy(snr_11bd,abs(mean(mse_11bd,2)),'Marker',markers(4))
semilogy(snr_11bd,abs(mean(mse_11bd_noW,2)),'Marker',markers(5))
semilogy(snr_UWcp,abs(mean(mse_UWcp,2)),'Marker',markers(6))
semilogy(snr_UWcp,abs(mean(mse_UWcp_noW,2)),'Marker',markers(7))
legend('UW with Wiener', 'UW w/o Wiener', '11p - preamble only', '11bd with Wiener', '11bd w/o Wiener','cpUW with Wiener', 'cpUW w/o Wiener')
ylabel('CE MSE')
xlabel('E_s/N_0') 

% SNR per data bit
EbNo_UW = snr_UW - 10*log10(mu*coderate);
EbNo_UWcp = snr_UWcp - 10*log10(mu*coderate);
EbNo_11p = snr_11p - 10*log10(mu*coderate);
EbNo_11bd = snr_11bd - 10*log10(mu*coderate);

figure; 
semilogy(EbNo_UW,abs(mean(ber_UW,2)),'Marker',markers(1))
hold on; 
semilogy(EbNo_UW,abs(mean(ber_UW_noW,2)),'Marker',markers(2))
semilogy(EbNo_11p,abs(mean(ber_CP,2)),'Marker',markers(3))
semilogy(EbNo_11bd,abs(mean(ber_11bd,2)),'Marker',markers(4))
semilogy(EbNo_11bd,abs(mean(ber_11bd_noW,2)),'Marker',markers(5))
semilogy(EbNo_UWcp,abs(mean(ber_UWcp,2)),'Marker',markers(6))
semilogy(EbNo_UWcp,abs(mean(ber_UWcp_noW,2)),'Marker',markers(7))
legend('UW with Wiener', 'UW w/o Wiener', '11p - preamble only', '11bd with Wiener', '11bd w/o Wiener','cpUW with Wiener', 'cpUW w/o Wiener')
ylabel('BER')
xlabel('E_b/N_0') 


spacesize_uw = (UW_length + N)*noblk;
spacesize_uwcp = (UW_length_withcp + N)*noblk;
spacesize_11p = (Ncp + N)*noblk;
spacesize_11bd = (Ncp+N)*(noblk+no_midambles);
BW_eff = length(data_subcar_indexes)/N;

spectral_eff_uw = BW_eff*(mu*noblk*N)/spacesize_uw;
spectral_eff_uwcp = BW_eff*(mu*noblk*N)/spacesize_uwcp;
spectral_eff_11p = BW_eff*(mu*noblk*N)/spacesize_11p;
spectral_eff_11bd = BW_eff*(mu*noblk*N)/spacesize_11bd;

fprintf('\nSpectral efficiencies for\n')
fprintf('UW: %4.2f bit/s/Hz \n',spectral_eff_uw)
fprintf('UW with CP: %4.2f bit/s/Hz \n',spectral_eff_uwcp)
fprintf('802.11p: %4.2f bit/s/Hz \n',spectral_eff_11p)
fprintf('802.11bd: %4.2f bit/s/Hz \n',spectral_eff_11bd)

if save_workspace
   save('last_run_workspace') 
   fprintf('\nSimulation data stored in \n%s\\last_run_workspace.mat\n',pwd)
end


