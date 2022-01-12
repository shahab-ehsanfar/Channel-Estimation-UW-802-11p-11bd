function [payload_withCP] = do_insert_CP(signal_inTime,Ncp)

payload_withCP = padarray(zeros(size(signal_inTime)),Ncp,0,'pre');
for blk=1:size(signal_inTime,3)
    for iT=1:size(signal_inTime,2)
        payload_withCP(:,1,blk) = [signal_inTime(end-Ncp+1:end,1,blk); signal_inTime(:,1,blk)];
    end
end

end

