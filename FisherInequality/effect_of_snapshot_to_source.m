% This script is used to check the effect of snr to MSE of source

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

testNumber = 20;
frequency=300e6;
c=3e8;
lambda=c/frequency;
spacing=0.5*lambda;
% signal to noise ratio, we define it by SNR = max(signal)/noise
snapshotSize=[1,5,10,25,50,100];
SNR = 20;
signal=[sqrt(2)*exp(1i*pi/4);sqrt(2)*exp(1i*pi/4);sqrt(2)*exp(1i*pi/4)]; % the signal 
% signal=1.5*exp(1i*1.45); % the signal 

signalVector=[real(signal);imag(signal)];

for snapIndex = 1 : length(snapshotSize)
    figure(snapIndex)
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
        antennaNumber = 12;

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
        % rng(1);
        [~,sigma] = my_awgn(ones(N_element,1)*yMaxElement,SNR);
        sigmaSquare = sigma^2; % the variance of the noise

        noise = my_awgn(ones(N_element,1)*yMaxElement,SNR)-ones(N_element,1)*yMaxElement; % complex numbers
        yWithNoiseComplex = yPure + noise;
        noiseVector = [real(noise);imag(noise)]; % real numbers form
        yWithNoise{i} = yPureVector + noiseVector;
        end

        initsigma2 = var(yWithNoise{1})*0.01;
        [time_l1_pair,recWeights,recIndex,sigma2,recSD,basis,selected,alpha,lambdas] = MultiTaskFastLaplacePair(entireMeasurementMatrix,yWithNoise,initsigma2,1e-6,lambda);

        color = [1 0 0; 0 1 0; 0 0 1; 0 1 1];
        recIndexCopy = recIndex;
        recIndexCopy(recIndex > N_discrete) = recIndex(recIndex > N_discrete) - N_discrete;
        pairAngle = unique(recIndexCopy);
        pairPowerSD = zeros(N_discrete,1);
        for i = 1 : length(pairAngle)
            if isempty(find(recIndex == pairAngle(i)))
                pairReal = 0;
            else
                pairReal = recSD(recIndex == pairAngle(i));
            end
            if isempty(find(recIndex == (pairAngle(i) + N_discrete)))
                pairImag = 0;
            else
                pairImag = recSD(recIndex == (pairAngle(i) + N_discrete));
            end
            pairPowerSD(pairAngle(i)) = pairReal^2 + pairImag^2;
        end
        recWeightsAvg = mean(recWeights,2);
        pairPower = recWeightsAvg(1:N_discrete).^2 + recWeightsAvg(N_discrete+1:end).^2;
        pairIdx = find(pairPower > 0);
    %     % plot error bar
    %     h = errorbar(pairIdx-91,pairPower(pairIdx),...
    %         pairPowerSD(pairIdx),'o','MarkerSize',5);
    %     h.CapSize = 10;
        % plot the mean
        hold on;
        h = plot(signalAngle,pairPower(signalAngle+91),'o','MarkerSize',5,'color','b');

        % the true signal power
        for i = 1 :length(signal)
            fh{i} = plot(signalAngle(i),abs(signal(i))^2,'*','LineWidth',1.5,'color',color(2,:)...
                ,'MarkerSize',10);
        end
        xlim([-95,95]);
        ylim([0,4]);
        xlabel('angle(theta)');
        ylabel('signal strength');
        lg = legend('reconstructed signal','true signal');
        set(gca,'fontsize',16)
        lg.FontSize = 14;
        % end of plotting resuls

        lowerBound = sigmaSquare/(2*antennaNumber*testNumber*0.1*snapshotNumber);
        actualWeight = zeros(size(recWeightsAvg,1),1);
        actualWeight(signalAngle+91) = 1;
        actualWeight(signalAngle+91+181) = 1;
        actualVariance = sum((actualWeight-recWeightsAvg).^2)/testNumber;
    end
    text(min(xlim)+5,max(ylim)-0.3,['lower bound = ',num2str(lowerBound)],'FontSize',11);
    text(min(xlim)+10,max(ylim)-0.5,['MSE = ',num2str(actualVariance)],'FontSize',11);
    fprintf('The lower bound for the variance : %f \n',lowerBound);
    fprintf('The MSE : %f \n',actualVariance);
end



