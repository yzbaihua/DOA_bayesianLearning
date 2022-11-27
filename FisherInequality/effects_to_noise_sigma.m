% This script is used to check the effect of snapshot to MSE of noise

% we used the bayesian compressive sensing method to recover the incident
% signal from the noised measurements

% Here we assume there are three clusters which are seperated by deltaDis(the
% default value is 100 lambda)

% The parameters you can modify are :
% 1) signal: the strength of the sources (complex column vector)
% 2) SNR: signal to noise ratio
% 3) signalAngle: the angles of the sources
% 4) snapshotNumber
clear variables;

testNumber = 2;
frequency=300e6;
c=3e8;
lambda=c/frequency;
spacing=0.5*lambda;
% signal to noise ratio, we define it by SNR = max(signal)/noise
snapshotSize=[10,50];
elementNumberSize=[6,9,12,15,18];
snrSize = [5,10,15,20,25,30];
signal=[sqrt(2)*exp(1i*pi/4);sqrt(2)*exp(1i*pi/4);sqrt(2)*exp(1i*pi/4)]; % the signal 
% signal=1.5*exp(1i*1.45); % the signal 

signalVector=[real(signal);imag(signal)];

for snapIndex = 1 : length(snapshotSize)
    for eleIndex = 1 : length(elementNumberSize)
        for snrIndex = 1 : length(snrSize)
            lowerBound = 0;
%             rng(1);
            for test = 1 : testNumber
                %     signalAngle=-1*randi([-90,90],1,2);%,-45,15,50]; % the angle of the signal
                signalAngle=[50,-20,15];%,-45,15,50]; % the an
                N_discrete=181; % the number of the possible angles
                possibleTheta=linspace(-90,90,N_discrete);
                
                % loc=([-8,-7,-3,-1.5,0,2,3.5,6,7,10]+8)*spacing; % the unit is lambda/2
                avgRecValueMean = zeros(2*181,1);
                avgRecValueSD = zeros(2*181,1);
                avgRecNoiseVariance = 0;
                snapshotNumber = snapshotSize(snapIndex);
                antennaNumber = elementNumberSize(eleIndex);
                
                for i = 1 : snapshotNumber
                    % loctaions contains three clusters, the clusters are seperated by 100
                    % lambda, and each cluster has three points, the locations are generated
                    % randomly
                    % rng('default')
                    %         % cluster configuration
                    %         deltaDis = 103*lambda;
                    %         loc =[3*rand(1,3), 3*rand(1,3)+deltaDis, 3*rand(1,3)+2*deltaDis,3*rand(1,3)+3*deltaDis];
                    
                    % general random array
                    deltaDis = 24*lambda;
                    for antennaIdx = 1 : antennaNumber
                        if antennaIdx == 1
                            loc = 6*lambda*rand;
                        else
                            loc = [loc,6*lambda*rand + (antennaIdx-1)*deltaDis];
                        end
                    end
                    
                    N_element=length(loc); % the number of the elements in the array
                    
                    % the measurement matrix contains the entrie potential angles
                    measurementMatrix = exp(-1i*loc'*2*pi*sind(possibleTheta)/lambda);
                    % covert it to the real number form
                    entireMeasurementMatrix{i} = ...
                        [real(measurementMatrix),-imag(measurementMatrix);imag(measurementMatrix),real(measurementMatrix)];
                    
                    
                    % we want to check the effect of position of elements to estimated results.
                    normObs = entireMeasurementMatrix{i};
                    prodObs = normObs'*normObs;
                    sumDiag = sum(diag(prodObs/norm(prodObs))); % !!! accuracy is proportional to that
                    
                    
                    % to get the actual received signal
                    actualMeasurementMatrix = exp(-1i*loc'*2*pi*sind(signalAngle)/lambda);
                    yPure = actualMeasurementMatrix * signal; % the received signal without noise, complex numbers
                    yPureVector = [real(yPure);imag(yPure)]; % real numbers form
                    yMaxElement = max(yPure); % the maxmium received signal is picked up for the SNR
                    
                    
                    % to control the noise
                    % rng('default');
                    
                    [~,sigma] = my_awgn(ones(N_element,1)*yMaxElement,snrSize(snrIndex));
                    sigmaSquare = sigma^2; % the variance of the noise
                    
                    noise = my_awgn(ones(N_element,1)*yMaxElement,snrSize(snrIndex))-ones(N_element,1)*yMaxElement; % complex numbers
                    yWithNoiseComplex = yPure + noise;
                    noiseVector = [real(noise);imag(noise)]; % real numbers form
                    yWithNoise{i} = yPureVector + noiseVector;
                end
                
                initsigma2 = var(yWithNoise{1})*0.01;
                [time_l1_pair,recWeights,recIndex,sigma2,recSD,basis,selected,alpha,lambdas] = MultiTaskFastLaplacePair(entireMeasurementMatrix,yWithNoise,initsigma2,1e-6,lambda);
                
                lowerBound = lowerBound + sigmaSquare^2/(2*antennaNumber^2*testNumber*snapshotNumber);
                actualSigmaSquareVariance = sum((sigma2-sigmaSquare).^2)/testNumber;
            end
            lowerBoundResults(snrIndex,(snapIndex-1)*length(elementNumberSize) + eleIndex) = lowerBound;
            MSEresults(snrIndex,(snapIndex-1)*length(elementNumberSize) + eleIndex) = actualSigmaSquareVariance;
            fprintf('L: %d, Element: %d, SNR: %d, The lower bound for the variance : %f \n',snapshotNumber, antennaNumber, snrSize(snrIndex), lowerBound);
            fprintf('The MSE : %f \n',actualSigmaSquareVariance);
        end
    end
end