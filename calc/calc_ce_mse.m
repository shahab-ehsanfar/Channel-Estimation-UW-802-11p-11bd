function [mse] = calc_ce_mse(config,h_genie,H_hat_out)
% Calculate the channel estimation MSE
%
% TU Chemnitz
% Shahab Ehsanfar 

N = config.N; % FFT size
nT = config.nT; % number of Tx antennas
nR = config.nR; % number of Rx antennas
noblk = config.noblk; % number of transmitted blocks/symbols. 

mse = zeros(noblk,1);
for blk = 1:noblk
    
    H_genie = zeros(N,nT,nR);
    for iR = 1:nR
        H_genie(:,:,iR) = fft(h_genie{blk}(:,:,iR),N);
    end
    
    mse(blk) = mean(abs(H_genie(:) - H_hat_out{blk}(:) ).^2);
    
end

mse = mean(mse);
end

