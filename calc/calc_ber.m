function [ber] = calc_ber(config,receivedBits)
% Calculate the Bit Error Rate

noblk = config.noblk;
b_blk = config.payload.b_blk;

ber=zeros(noblk,1);
for blk=1:noblk
   ber(blk) = mean( receivedBits{blk} ~= b_blk{blk});
end

ber= mean(ber);

end

