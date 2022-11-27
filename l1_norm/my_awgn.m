function [desiredSignal,noiseDeviation] = my_awgn( signal, SNR)
% signal: the actual signal without noise
% SNR : signal to noise ratio (dB)
% desiredSignal: the signal after the noise is added
% rng('default');
l = numel(signal);
snr_value = 10^(SNR/10); % convert to linear scale
signal_energy = sum(abs(signal(:)).^2)/l;
noise_energy = signal_energy/snr_value;
if(isreal(signal))
    noiseDeviation = sqrt(noise_energy);
    noiseSignal = noiseDeviation * randn(size(signal,1),size(signal,2));
else
    noiseDeviation = sqrt(noise_energy/2);
    noiseSignal = noiseDeviation * (randn(size(signal,1),size(signal,2))...
                  +1i*randn(size(signal,1),size(signal,2)));
end
desiredSignal = signal+noiseSignal;

end

