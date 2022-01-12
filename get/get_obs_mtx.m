function [config] = get_obs_mtx(config,UW_inTime_1,UW_length,numTaps,nT)
% This function has two sections
% 1) Creating the observation matrices for MMSE Channel estimation
% 2) Creating the observation matrices for MMSE equalization with imperfect channel knowledge
%
% TU Chemnitz
% Shahab Ehsanfar

%% Construct the observation matrix for channel estimation
N_preamble = UW_length*2;
% preamble
xp_iT = [UW_inTime_1; UW_inTime_1];

fft_path_delays_Np = 1:numTaps; % If the channel has non-zero taps at the end of the block, this needs to be updated. 
fft_path_delays_2Np = fft_path_delays_Np; % If the channel has non-zero taps at the end of the block, this needs to be updated. 

xp_NpreambleminusL_L_iT = zeros(N_preamble-numTaps+1,numTaps,nT);
xp_NpminusL_L_iT = zeros(UW_length-numTaps+1,numTaps,nT);
for iT = 1:nT
    % For secondary channel estimation based on UW
    xp_toeplitz = toeplitz(xp_iT(UW_length+1:end,iT),[xp_iT(1,iT) zeros(1,UW_length-1)]);
    xp_NpminusL_L_iT(:,:,iT) = xp_toeplitz(numTaps:end,fft_path_delays_Np);
    
    % For primary channel estimation based on preamble
    xp_toeplitz_preamble = toeplitz(xp_iT(:,iT),[xp_iT(1,iT) zeros(1,N_preamble-1)]);
    xp_NpreambleminusL_L_iT(:,:,iT) = xp_toeplitz_preamble(numTaps:end,fft_path_delays_2Np);
end
if nT == 4 % This part still needs to be automized 
   xp_NpminusL_L = [xp_NpminusL_L_iT(:,:,1) xp_NpminusL_L_iT(:,:,2) xp_NpminusL_L_iT(:,:,3) xp_NpminusL_L_iT(:,:,4)];
   xp_NpreambleminusL_L = [xp_NpreambleminusL_L_iT(:,:,1) xp_NpreambleminusL_L_iT(:,:,2) ...
                                xp_NpreambleminusL_L_iT(:,:,3) xp_NpreambleminusL_L_iT(:,:,4)];
elseif nT == 1
   xp_NpminusL_L = [xp_NpminusL_L_iT(:,:,1)];
   xp_NpreambleminusL_L = [xp_NpreambleminusL_L_iT(:,:,1)];
end
config.preamble.xp_NpminusL_L = xp_NpminusL_L;
config.preamble.xp_NpreambleminusL_L = xp_NpreambleminusL_L;
config.preamble.xp_iT = xp_iT;

%% Construct observation matrices for equalization with imperfect channel knowledge
% Generally the idea of these matrices come from GFDM where A is different
% than F^H. So feel free to reduce the complexity by removing extra
% identity matrices (which are specific for OFDM)

N = config.N;
fft_path_delays_Nd = 1:config.numPaths;

F_d = dftmtx(N)/sqrt(N);
A_ofdm = F_d';

F_L_d = F_d(:,fft_path_delays_Nd);
A_uw = A_ofdm; 
R_dd = eye(N); 
FA_Rdd_AF_tx1 = F_d*(A_uw*A_uw')*F_d';


config.matrices.A_uw = A_uw;
config.matrices.FA = F_d*A_uw;
config.matrices.F_L_d = F_L_d;
config.matrices.R_dd = kron(eye(nT),R_dd);
config.matrices.FA_Rdd_AF_tx1 = FA_Rdd_AF_tx1;


end

