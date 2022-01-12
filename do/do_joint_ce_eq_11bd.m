function [ d_hat, H_hat_out, Rx_matrices ] = do_joint_ce_eq_11bd(config, snr, y, h_hat_primary, R_hh_err_primary, noblk, wiener_flag)
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
Ncp = config.Ncp;
numTaps = config.numTaps; % number of channel taps
fd = config.fd; % maximum Doppler shift
Fs = config.Fs; % Sampling frequency
n_H = config.n_H_11bd; % sample index where CE is at its best
N_Nuw = config.N_Nuw; % N_Nuw = N (data) + Nuw = N (data) + N_preamble/2
midamble_inFreq = config.midamble_inFreq; 
midamble_spacing = config.midamble_spacing; % time domain midamble symbol spacing
no_midambles = config.no_midambles; % number of midambles
noblk = config.noblk;



processed_blks = 0; % initialization - number of equalized packets
midamble_indx = 0;

n_Hd = Ncp + N/2; 

H_hat_out = cell(no_midambles,1);
% h_hat_out{1} = h_hat_primary;
h_hat_all = h_hat_primary;

h_hat_in = h_hat_primary;
R_hh_err_in = R_hh_err_primary;

% separate the midamble symbols for channel estimation
y_midambles = y(:,:,midamble_spacing:midamble_spacing:end); 

d_hat = cell(noblk,1);
while processed_blks <= (noblk+no_midambles)          
    midamble_indx = midamble_indx+1;
    
    % Secondary Channel Estimation based on a single UW        
    if (midamble_indx <= no_midambles)
        h_hat_secondary = zeros(numTaps,nT,nR);
        R_hh_err_2nd = zeros(nT*numTaps,nT*numTaps,nR);
        for iR = 1:nR           
            h_hat_secondary(:,:,iR) = do_fft_channel_estimation(config,y_midambles(:,:,midamble_indx),midamble_inFreq(:,:,midamble_indx),snr,0,'low complexity approx. LMMSE - UW');
            R_hh_err_2nd(:,:,iR) = 1e-5*eye(nT*numTaps);
        end
    end
    

    % WIENER FILTERING
    if ( (wiener_flag == 1 && midamble_indx>1) || (wiener_flag == 2 && midamble_indx <= no_midambles))
        h_all_blks = zeros(numTaps,nT,nR,midamble_indx+midamble_spacing+1);
        Covariance = zeros(numTaps,numTaps,nT,nR,midamble_indx+midamble_spacing+1);
        for iR = 1:nR
            for iT = 1:nT
                h_hat_wiener_in = cat(2, squeeze(h_hat_all(:,iT,iR,:)), h_hat_secondary(:,iT,iR));
                R_hh_err_2nd_iT = R_hh_err_2nd((iT-1)*numTaps+1:iT*numTaps,(iT-1)*numTaps+1:iT*numTaps,iR);
                [h_all_blks(:,iT,iR,:), Covariance(:,:,iT,iR,:)] = do_wiener_smoodiction_11bd(h_hat_wiener_in, R_hh_err_2nd_iT, N+Ncp, midamble_spacing,n_H, n_Hd, Fs, fd, numTaps);
            end
        end
        h_hat_in = h_all_blks(:,:,:,end-midamble_spacing+1:end-1);
        
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
        
        
    elseif (wiener_flag == 2 && midamble_indx > no_midambles)
        h_hat_secondary = zeros(numTaps,nT,nR);
        h_all_blks = zeros(numTaps,nT,nR,midamble_indx+midamble_spacing);
        Covariance = zeros(numTaps,numTaps,nT,nR,midamble_indx+midamble_spacing);
        for iR = 1:nR
            for iT = 1:nT
                h_hat_wiener_in = cat(2, squeeze(h_hat_all(:,iT,iR,:)), h_hat_secondary(:,iT,iR));
                R_hh_err_2nd_iT = R_hh_err_2nd((iT-1)*numTaps+1:iT*numTaps,(iT-1)*numTaps+1:iT*numTaps,iR);
                [h_all_blks(:,iT,iR,:), Covariance(:,:,iT,iR,:)] = do_wiener_smoodiction_11bd(h_hat_wiener_in, R_hh_err_2nd_iT, N+Ncp, midamble_spacing,n_H, n_Hd, Fs, fd, numTaps);
            end
        end
        h_hat_in = h_all_blks(:,:,:,end-midamble_spacing+2:end);

        for iR = 1:nR
            if nT == 4
            R_hh_err_in(:,:,iR) = blkdiag((Covariance(:,:,1,iR,end)),(Covariance(:,:,2,iR,end)), ...
                (Covariance(:,:,3,iR,end)),(Covariance(:,:,4,iR,end)));
            elseif nT == 1
                R_hh_err_in(:,:,iR) = Covariance(:,:,1,iR,end);
            end
        end 
        
    elseif wiener_flag == 3 % Genie-Aided receiver with perfect CSI knowledge 
        h_hat_in = h_hat_primary{midamble_indx};
        R_hh_err_in = zeros(nT*numTaps, nT*numTaps, nR);
    else
%         if midamble_indx <= no_midambles
%             h_hat_in = repmat((h_hat_all(:,:,:,end) + h_hat_secondary)/2,[1 1 1 3]);
%             h_hat_all = cat(4, h_hat_all, h_hat_secondary);
%             R_hh_err_in = R_hh_err_2nd;
%         else
            h_hat_all = cat(4, h_hat_all, h_hat_secondary);
            h_hat_in = repmat(h_hat_all(:,:,:,end),[1 1 1 3]);            
%         end
        
        
    end
    
    
    for blk=1:(midamble_spacing-1)
        processed_blks = processed_blks +1;
        if processed_blks < noblk+no_midambles
            
            if mod(processed_blks,midamble_spacing)==0
               processed_blks = processed_blks + 1; 
            end
            
            yd = y(:,:,processed_blks);
            
            % FFT Operation
            Y_all = fft_u(yd);
            Y_all_v = Y_all(:);
            H_hat_in = zeros(N,nT,nR);
            for iR = 1:nR
                H_hat_in(:,:,iR) = fft(h_hat_in(:,:,iR,blk),N);
            end
            
            % Write out the channel frequency response
            for iR = 1:nR
                H_hat_out{processed_blks}(:,:,iR) = H_hat_in(:,:,iR);
            end            
            
            % MMSE Equalization
            [ d_hat{processed_blks} , ~, ~, R_dd_hat, ~] = do_mimo_mmse_equalization_cpofdm( config, snr, H_hat_in, R_hh_err_in, Y_all_v);
            Rx_matrices.R_dd_hat{processed_blks} = diag(R_dd_hat);          
            
            assert(length(find(abs(diag(R_dd_hat))<1e-8)) ~= length(diag(R_dd_hat)))
           
        end
        
    end

end

d_hat = d_hat(~cellfun('isempty',d_hat));
H_hat_out = H_hat_out(~cellfun('isempty',H_hat_out));
Rx_matrices.R_dd_hat = Rx_matrices.R_dd_hat(~cellfun('isempty',Rx_matrices.R_dd_hat));


end
