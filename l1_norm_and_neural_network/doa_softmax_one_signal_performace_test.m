clear variables;
frequency=300e6;
c=3e8;
lambda=c/frequency;
spacing=0.5*lambda;
% signal to noise ratio, we define it by SNR = max(signal)/noise
SNR=10; 

% signal=[1.1*exp(1i*1.05);0.8*exp(1i*0.3);1.5*exp(1i*1.3)]; % the signal 
signal=sqrt(2)*exp(1i*pi/4); % the signal 

% signalVector=[real(signal);imag(signal)];
deltaDis = 103*lambda;
loc =[1,2,3, [1,2,3]+deltaDis, [1,2,3]+2*deltaDis,[1,2,3]+3*deltaDis];
trainSize = 1810;
trainMatrixRealPart = zeros(length(loc),trainSize);
trainMatrixImagPart = zeros(length(loc),trainSize);
trainMatrixFull = zeros(2*length(loc),trainSize);
trainLabel = zeros(1,trainSize);
for trainIdx = 1 : trainSize
%     signalAngle=randi([-90,90],1,1);
    signalAngle= ceil((trainIdx-0.01)/10) - 91;
%     signal=sqrt(2*rand)*exp(1i*pi/(4*rand)); % the signal
    signal=sqrt(2)*exp(1i*pi/(4)); % the signal
    trainLabel(1,trainIdx) = signalAngle + 91;
    actualReceivedMatrix = exp(-1i*loc'*2*pi*sind(signalAngle)/lambda)*my_awgn(signal,SNR);
    trainMatrixRealPart(:,trainIdx) = real(actualReceivedMatrix);
    trainMatrixImagPart(:,trainIdx) = imag(actualReceivedMatrix); 
    trainMatrixFull(:,trainIdx) = [real(actualReceivedMatrix);imag(actualReceivedMatrix)];
end

inputSize = length(loc);
allLabels = unique(trainLabel);
numCls = length(allLabels);
softmaxLambda = 1e-4;
options.maxIter = 500;
% training part ---------
% predict real part
softmaxModelReal = softmaxTrain(inputSize,numCls,softmaxLambda,...
    trainMatrixRealPart,trainLabel,options);
saeSoftmaxOptThetaReal = softmaxModelReal.optTheta(:);
softmaxThetaReal = reshape(saeSoftmaxOptThetaReal(1:inputSize*numCls),numCls,inputSize);

% predict imag part
softmaxModelImag = softmaxTrain(inputSize,numCls,softmaxLambda,...
    trainMatrixImagPart,trainLabel,options);
saeSoftmaxOptThetaImag = softmaxModelImag.optTheta(:);
softmaxThetaImag = reshape(saeSoftmaxOptThetaImag(1:inputSize*numCls),numCls,inputSize);

% predict full part
softmaxModelFull = softmaxTrain(2*inputSize,numCls,softmaxLambda,...
    trainMatrixFull,trainLabel,options);
saeSoftmaxOptThetaFull = softmaxModelFull.optTheta(:);
softmaxThetaFull = reshape(saeSoftmaxOptThetaFull(1:2*inputSize*numCls),numCls,2*inputSize);
%----------------------------------------------


testNumber = 1;
correctNumber = 0;
thread = 0.9;
for i = 1 : testNumber
    signalAngleTest = unidrnd(180) - 90;
    signalTest = sqrt(2)*exp(1i*pi/4);
    noiseTest = my_awgn(signalTest,SNR)-signalTest;
    measurementTest = exp(-1i*loc'*2*pi*sind(signalAngleTest)/lambda)...
        *(signalTest+noiseTest);

    

    % predicting part ---------------
    measurementTestReal = real(measurementTest);
    measurementTestImag = imag(measurementTest);
    measurementTestFull = [measurementTestReal;measurementTestImag];
    [~, predReal] = max(softmaxThetaReal * measurementTestReal);
    [~, predImag] = max(softmaxThetaImag * measurementTestImag);
    [~, predFull] = max(softmaxThetaFull * measurementTestFull);
    
    [sortedProbability,sortIdx] = sort(exp(softmaxThetaFull * measurementTestFull)/sum(exp(softmaxThetaFull * measurementTestFull)),'descend');
    selectedIdx = 2;
    selectedAngle = sortIdx(1);
  
    while sum(sortedProbability(1:selectedIdx))/sum(sortedProbability) < thread 
        selectedAngle = [selectedAngle, sortIdx(selectedIdx)];
        selectedIdx = selectedIdx + 1;
    end
    if length(signalAngleTest) == 1
        if selectedAngle - 91 == signalAngleTest
            correctNumber = correctNumber + 1;
        end
    else
        for angleIdx = 1 : length(selectedAngle)
            if ismember(selectedAngle(angleIdx) - 91,signalAngleTest)
                correctNumber = correctNumber + 1;
            end
        end
    end
    angle = allLabels(predFull) - 91;
    
    % plot part 
%     figure()
%     plot(signalAngleTest, max(exp(softmaxThetaReal * measurementTestReal)),'*')
%     hold on
%     plot(allLabels - 91, exp(softmaxThetaReal * measurementTestReal))
%     figure()
%     plot(signalAngleTest, max(exp(softmaxThetaImag * measurementTestImag)),'*')
%     hold on
%     plot(allLabels - 91, exp(softmaxThetaImag * measurementTestImag))
% 
    figure()
    plot(signalAngleTest, max(exp(softmaxThetaFull * measurementTestFull))/sum(exp(softmaxThetaFull * measurementTestFull)),'*','MarkerSize',10)
    hold on
    plot(allLabels - 91, exp(softmaxThetaFull * measurementTestFull)/sum(exp(softmaxThetaFull * measurementTestFull)),'LineWidth',1.5);
    xlabel('Angles','FontSize',14)
    ylabel('Probability','FontSize',14)
end
accuracyPercentage = correctNumber/(testNumber*length(signalAngleTest));