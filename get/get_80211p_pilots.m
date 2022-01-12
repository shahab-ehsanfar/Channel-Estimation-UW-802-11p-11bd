function [pilots_inTime] = get_80211p_pilots(config)
% Gives out the pilots signal in time domain
% Pilots are PN sequences mapped to QPSK signals

N = config.N;
nT = config.nT; 
noblk = config.noblk;
pilot_subcar_indexes = config.pilot_subcar_indexes;

% PN sequence initialization
pnSequence = comm.PNSequence;
pnSequence.VariableSizeOutput = 1;
pnSequence.Polynomial = 'z^16 + z^13 + z + 1';
pnSequence.MaximumOutputSize = [2^16-1 1];
pnSequence.InitialConditions = [0 1 1 0 1 0 1 0 0 1 0 1 0 0 1 1];

% Sequence generation
PN_sequence = pnSequence(4*2*nT*noblk); % 2 bits for each pilot
PN_sequence_QPSK = qammod(PN_sequence,4,'InputType','bit','UnitAveragePower',true);

% Pilot allocation according to 802.11p subcarriers
pilots_inFreq = zeros(N,nT,noblk);
pilots_inTime = zeros(N,nT,noblk);
for blk= 1:noblk
    for iT = 1:nT
       pilots_inFreq(pilot_subcar_indexes,iT,blk) = PN_sequence_QPSK(((blk-1)*nT*4+(iT-1)*4+1:(blk-1)*nT*4+(iT-1)*4+4));
    end
    pilots_inTime(:,:,blk) = ifft_u(pilots_inFreq(:,:,blk));
end


end




