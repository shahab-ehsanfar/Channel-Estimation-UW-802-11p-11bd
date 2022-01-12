function [payload_withUW] = do_insert_UW(payload_inTime,UW_inTime)

if nargin==3
    mode='Circular Unique Words';
else
    mode='Identical Unique Words';
end

if strcmp(mode, 'Identical Unique Words')
    payload_withUW = padarray(zeros(size(payload_inTime)),length(UW_inTime),0,'pre');
    for blk=1:size(payload_inTime,3)
        for iT=1:size(payload_inTime,2)
            payload_withUW(:,1,blk) = [payload_inTime(:,1,blk); UW_inTime];
        end
    end
elseif strcmp(mode,'Circular Unique Words')
    % To be defined
end

end

