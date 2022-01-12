function [estimated_CFR,R_hh_err] = do_fft_channel_estimation(config,rx_signal,ref_signal,snr,freq_avrg_flag,method)
% Low complexity FFT based Channel Estimation for Wifi-based preamble
%
% snr (in dB)
%
% Shahab Ehsanfar
% TU Chemnitz

switch method
    case 'low complexity approx. LMMSE' % should be extended to MIMO
        T1_length = config.T1_length;
        noisVar = 10^(-0.1*snr);
        
        % Average the two copies of the preamble
        averaged_measurement = 0.5*(rx_signal(1:T1_length) + rx_signal(T1_length+1:end));
        
        % Estimate the channel in frequency domain
        rx_signal_inFreq = fft_u(averaged_measurement);
        raw_estimated_CFR = rx_signal_inFreq.*(conj(ref_signal)./(ref_signal.*conj(ref_signal) + noisVar));
        
        % Frequency domain averaging over 3 subcarriers
        if freq_avrg_flag
            filter_inFreq = rcosdesign(0.5,1,4);
            filter_inFreq(1:floor(end/2)) = [];
            filter_inFreq = filter_inFreq/sum(filter_inFreq);
            estimated_CFR = filter(filter_inFreq,1,[raw_estimated_CFR; zeros(2,1)]);
            estimated_CFR = estimated_CFR(2:end-1);
        else
            estimated_CFR = raw_estimated_CFR;
        end
        
        R_hh_err = [];
        
    case 'low complexity approx. LMMSE - UW' % should be extended to MIMO
        T1_length = config.T1_length;
        noisVar = 10^(-0.1*snr);
        numPaths = config.numPaths;
                
        % Estimate the channel in frequency domain
        rx_signal_inFreq = fft_u(rx_signal);
        raw_estimated_CFR = rx_signal_inFreq.*(conj(ref_signal)./(ref_signal.*conj(ref_signal) + noisVar));
        
        % Frequency domain averaging over 3 subcarriers
        if freq_avrg_flag
            filter_inFreq = rcosdesign(0.5,1,4);
            filter_inFreq(1:floor(end/2)) = [];
            filter_inFreq = filter_inFreq/sum(filter_inFreq);
            estimated_CFR = filter(filter_inFreq,1,[raw_estimated_CFR; zeros(2,1)]);
            estimated_CFR = estimated_CFR(2:end-1);
        else
            estimated_CFR = raw_estimated_CFR;
        end
        
        estimated_raw_CIR = ifft(estimated_CFR);
        estimated_CIR = estimated_raw_CIR(1:numPaths);
        
        estimated_CFR = estimated_CIR; 
        
        R_hh_err = [];
                    
    case 'exact LMMSE'
        T1_length = config.T1_length;
        R_hh = config.preamble.R_hh;
        nT = config.nT;
        numPaths = config.numPaths;
        
        % must be extended to MIMO
        X = diag(ref_signal(:,1));
        size1_X = size(X,1);
                        
        F_L_d = dftmtx(size1_X);
        F_L_d = F_L_d(:,1:numPaths);
        
        noisVar = 10^(-0.1*snr);
        
        % For primary channel estimation, average the T1 and T2
        if length(rx_signal) == 2*T1_length           
           rx_signal = 0.5*(rx_signal(1:T1_length) + rx_signal(T1_length+1:end));
        end
               
        MMSE_estimator = R_hh*F_L_d'*X'/( X*F_L_d*R_hh*F_L_d'*X' + noisVar*eye(size1_X) );
        estimated_CFR = MMSE_estimator*fft_u(rx_signal);
        
        R_hh_err = R_hh - MMSE_estimator*X*F_L_d*R_hh;
        
          

end

end

