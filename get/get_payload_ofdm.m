function [ payload, config] = get_payload_ofdm(config, noblk, mode)

N_data = config.N_data; % Number of data subcarriers
data_subcar_indexes = config.data_subcar_indexes;
nT = config.nT; % Number of Tx antennas
nR = config.nR; % Number of Rx antennas
mu = config.mu; % Modulation order
N = config.N; % FFT size
trellis= config.trellis; % trellis for Turbo encoder
interleaver_size = config.interleaver_size; % Size of code interleaver (before encoder)



if strcmp(mode,'dummy')
    payload = (randn(N_data, nT, noblk)+1i*randn(N_data, nT, noblk))/sqrt(2);    
    
elseif strcmp(mode,'coded')
    
    dd = cell(noblk,1);

    Dd = zeros(N,nT,noblk);
    xd_iT = zeros(N,nT,noblk);   
        
    b_blk = cell(noblk,1);
    dd_nT = cell(noblk,1);
    bc_nT_blk = cell(noblk,1);
    code_interleaver = cell(noblk,1);

    for blk = 1:noblk
        % Random interleavers
        code_interleaver{blk} = randperm(interleaver_size);
        
        % Random binary bits
        b_blk{blk} = round(rand([1 interleaver_size])).';
        
        % Channel Encoder
        turboEnc = comm.TurboEncoder(trellis,code_interleaver{blk});
        bc_nT_blk{blk} = turboEnc(b_blk{blk});
        
        % QAM Mapper
        dd_nT{blk} = qammod(bc_nT_blk{blk},2^mu,'InputType','bit','UnitAveragePower',true);
                
        % Spatial multiplexing in case of nT>=2
        dd{blk} = reshape(dd_nT{blk},[N_data nT]);
                
        for iT = 1:nT
            % Subcarrier allocation
            Dd(data_subcar_indexes,iT,blk) = dd{blk}(:,iT);
            
            % OFDM Modulation
            xd_iT(:,iT,blk) = ifft_u(Dd(:,iT,blk));
            
        end
    end
    config.matrices.Dd_indx = logical(abs(Dd(:,:,1)));
    
    [on_bins, on_bins1, off_bins, on_bins_single] = get_on_bins_ofdm(config, nT,nR);
    config.off_bins = off_bins;
    config.on_bins = on_bins;
    config.on_bins1 = on_bins1;
    config.on_bins_single = on_bins_single;
    
    config.payload.code_interleaver = code_interleaver;
    config.payload.b_blk = b_blk;
    config.payload.bc_nT_blk = bc_nT_blk;
    config.payload.dd_nT = dd_nT;   

    
    payload = xd_iT;
end



end

