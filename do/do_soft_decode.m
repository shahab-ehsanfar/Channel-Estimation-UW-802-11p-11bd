function [receivedBits_UW, d_hat_soft, Heq] = do_soft_decode(config, d_hat, Rx_matrices)
% Decode the data through soft-decision demodulation
%
% Shahab Ehsanfar
% TU Chemnitz

noblk = config.noblk;
code_interleaver = config.payload.code_interleaver;
mu = config.mu;
trellis = config.trellis;
data_subcar_indexes = config.data_subcar_indexes;

d_hat_soft = cell(noblk,1);
Heq = cell(noblk,1);
receivedBits_UW = cell(noblk,1);
for blk=1:noblk
    % Calculate the equivalent noise variance for each subcarrier
    R_dd_hat = Rx_matrices.R_dd_hat{blk};
    noisVar = (real(sqrt( (1./(R_dd_hat(data_subcar_indexes))) -1 )) );
    
    assert(isempty(find(isinf(noisVar),1)))
    
    % Soft-out demodulation
    demodSignal = qamdemod(d_hat{blk}(data_subcar_indexes),2^mu,'OutputType','llr','UnitAveragePower',true,'NoiseVariance',noisVar);
    
    % Turbo decoder
    turboDec = comm.TurboDecoder(trellis,code_interleaver{blk},'NumIterations',8);
    receivedBits_UW{blk} = turboDec(-demodSignal); % Demodulated signal is negated
end

end

