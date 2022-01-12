function [ x_iT ] = get_uw_sequence( UW_length, iT, type, PrimeNr)
% UW_length: length of unique word
% PrimeNr: The Prime number for length of Polyphase sequence P^2


if strcmp(type,'Zadoff-Chu')
    % Useful only for SISO
    % Not suitable for MIMO 
    x_iT = zadoffChuSeq(1,UW_length/2+1); 
    x_iT(end) = []; 
    x_iT = [x_iT; x_iT];
    
elseif strcmp(type,'Polyphase')
    
    P = PrimeNr;  
    assert(isprime(P),'Wrong input: P must be a prime number');
    
    x_iT_polyph = zeros(P^2,1);
    for n0 = 0:P-1
        for n1 = 0:P-1
            n = n0*P+n1;
            x_iT_polyph(n+1) = exp(1i*2*pi*iT*n0*n1/P);
        end        
    end
    
    % Virtual subcarrier insertion
    x_iT_single = ifft_u([zeros(floor((UW_length - P^2)/2),1); fft_u(x_iT_polyph); zeros(ceil((UW_length - P^2)/2),1)]);
    
    x_iT = x_iT_single*sqrt(UW_length/P^2);
  
end


end

