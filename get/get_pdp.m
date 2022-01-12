% Copyright (c) 2020 TU Chemnitz
% All rights reserved.
% See accompanying license.txt for details.
%

function [P] = get_pdp(nT,nR,numPaths,AveGaindb,profile)
% Get the Power-Delay-Profile
%  
%  @param nT: number of Tx antennas
%  @param nR: number of Rx antennas
%  @param numPaths: Number of channel paths
%  @param AveGain: Average channel gains
%  @param profile: string 'exp decaying' or 'uniform'

switch profile
    case 'exp decaying'
        P = zeros(numPaths,nT,nR);
        for iR = 1:nR
            for iT = 1:nT
                AveGain = 10.^(AveGaindb/10);
                P(:,iT,iR) = AveGain./sum(AveGain); %a./sqrt(sum(a.^2));
            end
        end
    case 'uniform'
        P = zeros(numPaths,nT,nR);
        for iR = 1:nR
            for iT = 1:nT
                AveGain = 1/numPaths;
                P(:,iT,iR) = AveGain./sum(AveGain); %a./sqrt(sum(a.^2));
            end
        end
end
end

