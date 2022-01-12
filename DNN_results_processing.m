clc; clear variables; 

print_logs = 1;

 algo = 'cpUW_noW'; h_genie_name = 'h_genie_UWcp_dnn';
% algo = 'i11bd_noW'; h_genie_name = 'h_genie_CP_11bd_dnn';
%hlsizes = '453045'; % hl1 hl2 hl3 ...
hlsizes = '351835'; % hl1 hl2 hl3 ...

load('Raw_DNN_dataset_fd_1600_snr_1');
% load('Raw_DNN_dataset_snr_1');
load([algo '_DNN_Dataset_fd_1600_snr_1']);

i = 1; % Trained SNR index
testing_indices = DNN_Datasets.testing_indices;

N_Test_Frames = length(testing_indices);
% nUSC = DNN_Datasets.nUSC;
% nSym = DNN_Datasets.nSym;
nUSC = 64;
nSym = 120;
noblk = nSym; 

its_test = 300; % Number of simulation iterations dedicated to dNN test (~0.2 its)

PrimeNr = 5;
Ncp = 16;
UW_length_withcp = Ncp + PrimeNr^2;

numPaths = 8;


ber_UWcp_noW = zeros(length(snr),its_test);

for si = 1:length(snr)
    if si == 8
        [];
    end
    load(['Raw_DNN_dataset_fd_1600_snr_' int2str(si)]);
    load([algo '_DNN_' hlsizes '_Results_' int2str(si) '.mat']);    
    
    if strcmp(algo, 'cpUW_noW')
        y = y_uwcp;
    elseif strcmp(algo, 'i11bd_noW')
        y = y_cpofdm_11bd;
    end
    
    % Observations
    UW_noW_Testing_X = eval([algo '_DNN_' hlsizes '_test_x_',num2str(si)]);
    UW_noW_Testing_X = reshape(UW_noW_Testing_X(1:nUSC,:) + 1i*UW_noW_Testing_X(nUSC+1:2*nUSC,:), nUSC, nSym, N_Test_Frames);
    
    % Target vectors
    UW_noW_Testing_Y = eval([algo '_DNN_' hlsizes '_test_y_',num2str(si)]);
    UW_noW_Testing_Y = reshape(UW_noW_Testing_Y(1:nUSC,:) + 1i*UW_noW_Testing_Y(nUSC+1:2*nUSC,:), nUSC, nSym, N_Test_Frames);
    
    % Predictions
    UW_noW_DNN_Y = eval([algo '_DNN_' hlsizes '_corrected_y_',num2str(si)]);
    UW_noW_DNN_Y = reshape(UW_noW_DNN_Y(1:nUSC,:) + 1i*UW_noW_DNN_Y(nUSC+1:2*nUSC,:), nUSC, nSym, N_Test_Frames);
    
    H_hat_UW_dnn = cell([noblk 1]);
    h_genie_UW = cell([noblk 1]);
    H_hat_out_UW_noW = cell([noblk 1]);
    for i = 1:its_test
        if print_logs
            if mod(i,round(its/4)) == 0
                fprintf('at %s, running at iteration %d.\n',datestr(now,'HH:MM:SS'),i)
            elseif i == 1
                fprintf('Simulating SNR %d dB at max Doppler shift %d Hz. \nTime: %s\n',snr(si),config{i}.fd,datestr(now,'HH:MM:SS'))
            end
        end
               
        i_in = find(testing_indices==i);
        
        % Double check whether the target vectors are the original genie
        % channels
        h_test = ifft(UW_noW_Testing_Y(:,1,i_in));
        assert(mean(abs(h_test(1:8) - eval([h_genie_name '{1,i}'])))^2 < 1e-20);
        
        if strcmp(algo,'cpUW_noW')
            % fit the format to use the calc ce mse function
            for blk = 1:noblk
                H_hat_UW_dnn{blk} = UW_noW_DNN_Y(:,blk,i_in);          %  Result of dNN estimation
                h_genie_UW{blk} = ifft(UW_noW_Testing_Y(:,blk,i_in));
                %h_genie_UW{blk} = h_genie_UWcp_dnn{blk,testing_indices(i)};
                
                H_hat_out_UW_noW{blk} = UW_noW_Testing_X(:,blk,i_in);
                
                
                % CP Emulation at the receiver
                h_hat_in = ifft(H_hat_UW_dnn{blk});
                yd = do_emulate_cp(config{i},h_hat_in(1:numPaths),y{i},blk,'for_data_cpUW');
                
                % FFT Operation
                Y_all = fft_u(yd);
                Y_all_v = Y_all(:);
                H_hat_in = zeros(N,1,1);
                for iR = 1:1
                    H_hat_in(:,:,iR) = fft(h_hat_in(:,:,iR),N);
                end
                R_hh_err_in = 1e-6*eye(numPaths);
                
                % UW MMSE Equalization
                [d_hat{blk}, ~, ~, R_dd_hat, ~] = do_mimo_uw_mmse_equalization(config{i}, snr(si), H_hat_in, R_hh_err_in, Y_all_v);
                
                % Rx_matrices.R_dd_tilde{blk} = R_dd_tilde;
                
                Rx_matrices.R_dd_hat{blk} = diag(R_dd_hat);
                
                
            end
            [receivedBits, ~, ~] = do_soft_decode(config{i}, d_hat, Rx_matrices);
            clear Rx_matrices d_hat
            
        elseif strcmp(algo,'i11bd_noW')
            
            processed_blks = 0;
            midamble_indx = 0;
            no_midambles = config{1}.no_midambles;
            midamble_spacing = config{1}.midamble_spacing;
            blk = 0;
            
            while processed_blks <= (noblk+no_midambles)
                midamble_indx = midamble_indx+1;                
                
                for inter_midamble_blk=1:(midamble_spacing-1)                    
                    blk = blk + 1;
                    processed_blks = processed_blks +1;
                    if processed_blks < noblk+no_midambles
                        
                        if mod(processed_blks,midamble_spacing)==0
                            processed_blks = processed_blks + 1;
                        end
                        
                        H_hat_UW_dnn{blk} = UW_noW_DNN_Y(:,blk,i_in);          %  Result of dNN estimation
%                         H_hat_UW_dnn{blk} = UW_noW_Testing_X(:,blk,i_in); % debug
                        h_genie_UW{blk} = ifft(UW_noW_Testing_Y(:,blk,i_in));
                        %h_genie_UW{blk} = h_genie_UWcp_dnn{blk,testing_indices(i)};
                        
                        H_hat_out_UW_noW{blk} = UW_noW_Testing_X(:,blk,i_in);
                        
                        h_hat_in = ifft(H_hat_UW_dnn{blk});
                        
                        yd = y{i}(:,:,processed_blks);
                        
                        % FFT Operation
                        Y_all = fft_u(yd);
                        Y_all_v = Y_all(:);
                        H_hat_in = zeros(N,1,1);
                        for iR = 1:1
                            H_hat_in(:,:,iR) = fft(h_hat_in(1:numPaths),N);
                        end
                        R_hh_err_in = 1e-6*eye(numPaths);
                        
                        % MMSE Equalization
                        [ d_hat{blk} , ~, ~, R_dd_hat, ~] = do_mimo_mmse_equalization_cpofdm( config{i}, snr(si), H_hat_in, R_hh_err_in, Y_all_v);
                        Rx_matrices.R_dd_hat{processed_blks} = diag(R_dd_hat);
                        
                        assert(length(find(abs(diag(R_dd_hat))<1e-8)) ~= length(diag(R_dd_hat)))
                        
                    end
                end
            end
            d_hat = d_hat(~cellfun('isempty',d_hat));
            Rx_matrices.R_dd_hat = Rx_matrices.R_dd_hat(~cellfun('isempty',Rx_matrices.R_dd_hat));
            [receivedBits, ~, ~] = do_soft_decode(config{i}, d_hat, Rx_matrices);
            clear Rx_matrices d_hat
        end
        
        
        
        ber_UWcp_noW(si,i) = calc_ber(config{i},receivedBits);
        
        mse_X(si,i) = calc_ce_mse(config{i},h_genie_UW,H_hat_out_UW_noW);
        mse_dNN(si,i) = calc_ce_mse(config{i},h_genie_UW,H_hat_UW_dnn);
                        
    end
    if print_logs
        fprintf('Simulation of SNR %d ended.\n',snr(si))
        if si==length(snr)
            fprintf('Time: %s\n',datestr(now,'HH:MM:SS'))
        end
    end

end

snr_UWcp = snr +(UW_length_withcp + N)/N;

if length(snr)==1
    mean(mse_X)
    mean(mse_dNN)
else
    figure;
    semilogy(snr_UWcp,abs(mean(mse_X,2)))
    hold on
    semilogy(snr_UWcp,abs(mean(mse_dNN,2)))
end
    
figure; plot(abs(H_hat_UW_dnn{100}))
hold on; plot(abs(H_hat_out_UW_noW{100}))
hold on; plot(abs(UW_noW_Testing_Y(:,100,i_in)))
hold on; plot(abs(fft(eval([h_genie_name '{100,(i)}']),N)),'--')


EbNo_UWcp = snr_UWcp - 10*log10(config{1}.mu*config{1}.coderate);
figure
semilogy(EbNo_UWcp,abs(mean(ber_UWcp_noW,2)))

save('DNN_results_workspace.mat')
