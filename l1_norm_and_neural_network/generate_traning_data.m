% This script is used to generate the training data for certain locations
% input: 1) the locations of elements, 2) the angles (-90~90) 
% output: 1) the received signals, 2) the labels (i.e., the angles)

clear variables;
frequency=300e6;
c=3e8;
lambda=c/frequency;
% signal to noise ratio, we define it by SNR = max(signal)/noise
SNR=30; 
sampleResolution = 181; % the resolution of the training data
sampleNumber = 10; % the number of the samples in each angle
signalAngle = linspace(-90,90,sampleResolution);
signal = sqrt(2)*exp(1i*pi/4); % the signal 

signalVector=[real(signal);imag(signal)];
loc = [5.92966283554869,29.5769068345537,50.4570896006329,72.0020487733088,99.2452699365783,121.246383731976,145.315701529095,169.954837576405,192.575696499646,220.485201970430,244.491053744360,267.259796598511];
N_discrete=181; % the number of the possible angles
possibleTheta=linspace(-90,90,N_discrete);

yWithNoise = zeros(sampleResolution*sampleNumber,length(loc)*2);
yLabels = zeros(sampleResolution*sampleNumber,sampleResolution);
for i = 1 : length(signalAngle)
    incomingAngle = signalAngle(i);
    for j = 1 : sampleNumber
        % to get the actual received signal
        actualMeasurementMatrix = exp(-1i*loc'*2*pi*sind(incomingAngle)/lambda);
        yPure = actualMeasurementMatrix * signal; % the received signal without noise, complex numbers
        yPureVector = [real(yPure);imag(yPure)]; % real numbers form
        yMaxElement = max(yPure); % the maxmium received signal is picked up for the SNR



        [~,sigma] = my_awgn(ones(length(loc),1)*yMaxElement,SNR);
        sigmaSquare = sigma^2; % the variance of the noise

        noise = my_awgn(ones(length(loc),1)*yMaxElement,SNR)-ones(length(loc),1)*yMaxElement; % complex numbers
        yWithNoiseComplex = yPure + noise;
        noiseVector = [real(noise);imag(noise)]; % real numbers form
        yWithNoise((i-1)*sampleNumber+j,:) = yPureVector + noiseVector;
        yLabels((i-1)*sampleNumber+j,i) = 1;
    end
end

csvwrite('oneSourceData.csv',yWithNoise)
csvwrite('oneSourceLabelData.csv',yLabels)





