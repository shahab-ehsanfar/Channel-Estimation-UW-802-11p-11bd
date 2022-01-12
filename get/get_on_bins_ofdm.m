function [on_bins, on_bins1, off_bins, on_bins_single] = get_on_bins_ofdm(config, nT,nR)

N= config.N;

% Set of active data subcarriers for a single OFDM symbol 
on_bins_single = config.data_subcar_indexes;

% Set of active data subcarriers for all Tx antennas
on_bins1 = [];
for iT = 0:nT-1
   on_bins1 = [on_bins1; (on_bins_single+iT*N)]; 
end

% Set of active data subcarriers for all antennas
on_bins = [];
for iR =  0:nR-1
    on_bins = [on_bins; (on_bins1+iR*nT*N)];
end

% Set of inactive subcarriers
off_bins = 1:N*nR;
off_bins(on_bins1) = [];


end
