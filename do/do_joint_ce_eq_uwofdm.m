function [ d_hat, H_hat_out, Rx_matrices ] = do_joint_ce_eq_uwofdm(config, snr, y, h_hat_primary, R_hh_err_primary, noblk, wiener_flag)
% Joint Channel Estimation and Equalization for Unique Words
%
% Shahab Ehsanfar, TU Dresden
%
% wiener_flag = 0: Without wiener filtering
% wiener_flag = 1: Normal wiener filtering (smoothing and interpolation)
% wiener_flag = 2: (Reduces the performance, currently commented out) Wiener filtering with prediction of channel for future blocks (smoothing, interpolation and prediction)


nT = config.nT;
nR = config.nR;
N = config.N; 
UW_length = config.UW_length;
numTaps = config.numTaps; % number of channel taps
fd = config.fd; % maximum Doppler shift
Fs = config.Fs; % Sampling frequency
n_H = config.n_H; % sample index where CE is at its best
N_Nuw = config.N_Nuw; % N_Nuw = N (data) + Nuw = N (data) + N_preamble/2


B = 8; % Number of blocks to be processed for prediction (mainly for wiener_flag = 2)

% x_iT = zadoffChuSeq(1,UW_length/2+1); 
% x_iT(end) = []; 

n_Hd = UW_length + N/2; 

H_hat_out = cell(noblk,1);
% h_hat_out{1} = h_hat_primary;
h_hat_all = h_hat_primary;

h_hat_in = h_hat_primary;
R_hh_err_in = R_hh_err_primary;


d_hat = cell(noblk,1);
for blk = 1:noblk  
    
    
    yp = y(N+1:end,:,blk);
    
    % Secondary Channel Estimation based on a single UW
    h_hat_secondary = zeros(numTaps,nT,nR);
    R_hh_err_2nd = zeros(nT*numTaps,nT*numTaps,nR);
    if (wiener_flag ~= 2 || blk <= B)
        for iR = 1:nR           
            [h_hat_secondary(:,:,iR), R_hh_err_2nd(:,:,iR)] = do_mimo_uw_channel_estimation( config, snr, yp(:,iR), 'secondary_lmmse_skipL');
%             [h_hat_secondary(:,:,iR), R_hh_err_2nd(:,:,iR)] = do_fft_channel_estimation(config,yp(UW_length/2+1:end,iR),x_iT,snr,'true','exact LMMSE');
        end
    elseif blk > B
%         for iR = 1:nR
%             [h_hat_secondary(:,:,iR), R_hh_err_2nd(:,:,iR)] = do_mimo_uw_channel_estimation( config, snr, yp(:,iR), 'secondary_slmmse_skipL',h_hat_in(:,:,iR),R_hh_err_in(:,:,iR));
%         end
    end
    

    % WIENER FILTERING
    if ( (wiener_flag == 1 && blk>1) || (wiener_flag == 2 && blk < B))
        h_all_blks = zeros(numTaps,nT,nR,blk+2); % +2: (new secondary estimate) + (center of data symbol index)
        Covariance = zeros(numTaps,numTaps,nT,nR,blk+2);
        for iR = 1:nR
            for iT = 1:nT
                h_hat_wiener_in = cat(2, squeeze(h_hat_all(:,iT,iR,:)), h_hat_secondary(:,iT,iR));
                R_hh_err_2nd_iT = R_hh_err_2nd((iT-1)*numTaps+1:iT*numTaps,(iT-1)*numTaps+1:iT*numTaps,iR);
                [h_all_blks(:,iT,iR,:), Covariance(:,:,iT,iR,:)] = do_wiener_smoodiction_uw(h_hat_wiener_in, R_hh_err_2nd_iT, size(h_hat_wiener_in,2), N_Nuw, n_H, n_Hd, Fs, fd, numTaps);
            end
        end
        h_hat_in = h_all_blks(:,:,:,end-1);
        
        h_hat_all = cat(4, h_hat_all, h_hat_secondary);
        %     h_hat_in = h_hat_secondary;
        for iR = 1:nR
            if nT == 4
            R_hh_err_in(:,:,iR) = blkdiag((Covariance(:,:,1,iR,end)),(Covariance(:,:,2,iR,end)), ...
                (Covariance(:,:,3,iR,end)),(Covariance(:,:,4,iR,end)));
            elseif nT == 1
                R_hh_err_in(:,:,iR) = Covariance(:,:,1,iR,end);
            end
        end
        
        % Write out the channel frequency response (for blocks >= 2)
        for iR = 1:nR
            H_hat_out{blk}(:,:,iR) = fft(h_hat_in(:,:,iR),N); 
        end
    elseif (wiener_flag == 2 && blk >= B)
%         h_all_blks = zeros(numPaths,nT,nR,blk+3);
%         Covariance = zeros(numPaths,numPaths,nT,nR,blk+3);
%         for iR = 1:nR
%             for iT = 1:nT
%                 h_hat_wiener_in = cat(2, squeeze(h_hat_all(:,iT,iR,:)), h_hat_secondary(:,iT,iR));
%                 R_hh_err_2nd_iT = R_hh_err_2nd((iT-1)*numPaths+1:iT*numPaths,(iT-1)*numPaths+1:iT*numPaths,iR);
%                 [h_all_blks(:,iT,iR,:), Covariance(:,:,iT,iR,:)] = do_wiener_smoodiction_uw(h_hat_wiener_in, R_hh_err_2nd_iT, size(h_hat_wiener_in,2)+1, N, n_H, Fs, fd, numPaths);
%             end
%         end
%         h_hat_in = h_all_blks(:,:,:,end);
%         
%         h_hat_all = cat(4, h_hat_all, h_all_blks(:,:,:,end-1));
%         %     h_hat_in = h_hat_secondary;
%         for iR = 1:nR
%             R_hh_err_in(:,:,iR) = blkdiag((Covariance(:,:,1,iR,end)),(Covariance(:,:,2,iR,end)), ...
%                 (Covariance(:,:,3,iR,end)),(Covariance(:,:,4,iR,end)));
%         end
%         
%         h_hat_out{blk+1} = h_all_blks(:,:,:,end-1);
    elseif wiener_flag == 3 % Genie-Aided receiver with perfect CSI knowledge 
        h_hat_in = h_hat_primary{blk};
        R_hh_err_in = zeros(nT*numTaps, nT*numTaps, nR);
    else
        h_hat_in = (h_hat_all(:,:,:,end) + h_hat_secondary)/2;
        h_hat_all = cat(4, h_hat_all, h_hat_secondary);        
        R_hh_err_in = R_hh_err_2nd;
        
        % Write out the channel frequency response
        for iR = 1:nR
            H_hat_out{blk}(:,:,iR) = fft(h_hat_in(:,:,iR),N); 
        end
    end
    
    
    % CP Emulation at the receiver
    yd = do_emulate_cp(config,h_hat_in,y,blk,'for_data');
    
    % FFT Operation
    Y_all = fft_u(yd);
    Y_all_v = Y_all(:);
    H_hat_in = zeros(N,nT,nR);
    for iR = 1:nR
        H_hat_in(:,:,iR) = fft(h_hat_in(:,:,iR),N);
    end
    
    
    % UW MMSE Equalization
    [d_hat{blk}, ~, ~, R_dd_hat, ~] = do_mimo_uw_mmse_equalization(config, snr, H_hat_in, R_hh_err_in, Y_all_v);

        

        
    
    
   % Rx_matrices.R_dd_tilde{blk} = R_dd_tilde;
    
    Rx_matrices.R_dd_hat{blk} = diag(R_dd_hat);
    

end


end
