function [psd_filtered,f] = calc_psd(input_signal,NTotCarriers)
% Calculate the power spectral density of the input signal
input_signal0 = input_signal;


f = linspace(-NTotCarriers/2, NTotCarriers/2, 2*length(input_signal)+1); 
f = f(1:end-1)';
psd_exact = ((fftshift(abs(fft(input_signal, 2*length(input_signal))))).'/2).';
L = 100;
psd_filtered = abs(ifft(fft(psd_exact).*fft(ones(1,L),length(psd_exact)).')/sqrt(L));
plot_set = 1:50:length(f);
figure;
plot(f(plot_set), -59.0+mag2db(psd_filtered(plot_set)), 'b');
end

