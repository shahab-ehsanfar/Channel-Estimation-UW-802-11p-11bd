% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%


function [xch,h_genie] = do_mimo_channel(mode, x_mimo, h, numPaths, nT, nR, signal_length, noblk, data_middle_points)
% Apply channel convolution and noise to a signal X
%
% x_mimo: transmit signal
% h: matrix of channel imp. resp.
% nT: number of transmit antennas
% noblk: number of transmitted gfdm data blocks
% signal_length: including cp/uw length
%

switch mode
    case 'Filter'
        xch = zeros(signal_length,nR);
        xch_iRiT = zeros(signal_length,nT);        
        for iR = 1:nR
            for iT = 1:nT
                start_pnt = (iT-1)*numPaths+1;
                end_pnt = start_pnt + numPaths-1;
                                
                h_ii = h(start_pnt:end_pnt,iR);
                
                % Filter the signal with channel imp. resp.
                xch_iRiT(:,iT) = filter(h_ii,1,x_mimo(:,iT));
            end
            xch(:,iR) = sum(xch_iRiT,2);
        end
        
    case 'Filter time variant'
        xch = zeros(signal_length,nR);
        xch_iRiT = zeros(signal_length,nT);
        pathgains = cell(nT,nR);
        for iR = 1:nR
            for iT = 1:nT                                
                % Filter the signal with channel imp. resp.
                [xch_iRiT(:,iT),pathgains{iT,iR}] = h{iT,iR}(x_mimo(:,iT));
            end
            xch(:,iR) = sum(xch_iRiT,2);
        end
        
    case 'Toeplitz' % Slower because of matrix multiplication
        
        h_matrix = [];
        for iR = 1:nR
            h_matrix_iR = [];
            for iT = 1:nT
                start_pnt = (iT-1)*numPaths+1;
                end_pnt = start_pnt + numPaths-1;
            
                h_long = ifft(fft(h(start_pnt:end_pnt,iR),signal_length));
                
                % Create lower triangular channel matrix in time
                h_matrix_iRiT = toeplitz(h_long,[h_long(1); zeros(length(h_long)-1,1)]);                 
                h_matrix_iRiT(abs(h_matrix_iRiT) < 10e-16) = 0;
                
                % construct the mimo channel matrix in time
                h_matrix_iR = [h_matrix_iR h_matrix_iRiT];
            end
            h_matrix = [h_matrix; h_matrix_iR];
        end        
        xch = h_matrix*x_mimo(:);
    
end


% Extract and reorganize the Genie-aided channel knowledge
h_genie = cell(noblk,1);
for blk = 1:noblk
    h_genie{blk} = zeros(size(pathgains{1,1},2),nT,nR);
    for iR = 1:nR
        for iT = 1:nT
            h_genie{blk}(:,iT,iR) = pathgains{iT,iR}(data_middle_points(blk),:).';
        end
    end
end



end


















        
        
        
       
        