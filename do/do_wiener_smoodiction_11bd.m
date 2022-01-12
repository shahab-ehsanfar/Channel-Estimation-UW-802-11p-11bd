function [ h_all_blks, Covariance] = do_wiener_smoodiction_11bd(h_hat, R_hh_err, N_tot, midamble_spacing,n_H, n_Hd, Fs, fd, numPaths)
% Wiener filtering for joint Smoothing & Prediction for time-variant
% channels over multiple blocks
%
% Shahab Ehsanfar, TU Dresden
% n_H: sample index where channel estimation is at its best

% Primary channel estimation should have a double weight
H_in_time = [repmat(h_hat(:,1),[1 2]) h_hat(:,2:end)];

Ncp_long = 64;
N = 64;
preamble_blks = (Ncp_long + N/2):(Ncp_long+N):2*(Ncp_long+N);

midmbl_blk = [1 2:midamble_spacing:(size(H_in_time,2)-1)*midamble_spacing];
all_blk = [midmbl_blk(1:end-2) midmbl_blk(end-1):midmbl_blk(end)];

if mean(abs(H_in_time(:,end)))<1e-10 % If the last column is all zero, do prediction/extrapolation
    midmbl_blk(end)=[];
    all_blk(end)=[];
    H_in_time(:,end) = [];
end

[~, midmbl_indx] = intersect(all_blk,midmbl_blk,'stable'); 

R_f_indx = (midmbl_blk-1)*N_tot+ n_H; % sample index at midamble symbols
R_f_indx_all = (all_blk-1)*N_tot+ n_Hd; % sample index at all symbols

% R_f_indx_all = [R_f_indx_all(1:end-1) (noblk-2)*N_tot+n_Hd R_f_indx_all(end)];
% pilot_blk = [1:no_p_blk-1 no_p_blk+1];

% Channel time autocorrelation at pilot blocks
R_hp = zeros(length(midmbl_blk));
for i = 1:length(midmbl_blk)
    for j = 1:length(midmbl_blk)
        R_hp(i,j) = besselj(0,2*pi*( R_f_indx(i) - R_f_indx(j) )/Fs*fd);
    end
end

% Channel time autocorrelation at all blocks
no_wblk = length(all_blk);
R_h = zeros(no_wblk);
for i =1:no_wblk
    for j = 1:no_wblk
        R_h(i,j) = besselj(0,2*pi*( R_f_indx_all(i) - R_f_indx_all(j) )/Fs*fd);
    end
end

H_smoodicted_time = zeros(numPaths, no_wblk);
Covariance_l = zeros(no_wblk,no_wblk,numPaths);
Covariance = zeros(numPaths,numPaths,no_wblk);
for l=1:numPaths
   Wiener_mtx = (R_h(midmbl_indx,:)'/(R_hp+ R_hh_err(l,l)*eye(length(midmbl_blk))));
   H_smoodicted_time(l,:) = Wiener_mtx*H_in_time(l,:).';
   Covariance_l(:,:,l) = R_h - Wiener_mtx*R_h(midmbl_indx,:);
   % mse(l) = abs(mean(diag(Covariance(:,:,l))));
   for blk = 1:no_wblk
       Covariance(l,l,blk) = Covariance_l(blk,blk,l);
   end   
end
% figure; plot(mse)

h_all_blks = H_smoodicted_time; 



%Covariance = mean(Covariance_l,3);
 
end

