clear variables;
frequency=300e6;
c=3e8;
lambda=c/frequency;
spacing=0.5*lambda;
% signal to noise ratio, we define it by SNR = max(signal)/noise
SNR=30; 

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

signalAngleTest = [38, -35];
signalTest = [sqrt(2)*exp(1i*pi/4);sqrt(2)*exp(1i*pi/4)];
noiseTest = my_awgn(signalTest,SNR)-signalTest;
measurementTest = exp(-1i*loc'*2*pi*sind(signalAngleTest)/lambda)...
    *(signalTest+noiseTest);

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

% predicting part ---------------
measurementTestReal = real(measurementTest);
measurementTestImag = imag(measurementTest);
measurementTestFull = [measurementTestReal;measurementTestImag];
[~, predReal] = max(softmaxThetaReal * measurementTestReal);
[~, predImag] = max(softmaxThetaImag * measurementTestImag);
[~, predFull] = max(softmaxThetaFull * measurementTestFull);

angle = allLabels(predReal) - 91;
figure()
plot(signalAngleTest, max(exp(softmaxThetaReal * measurementTestReal)),'*')
hold on
plot(allLabels - 91, exp(softmaxThetaReal * measurementTestReal))
figure()
plot(signalAngleTest, max(exp(softmaxThetaImag * measurementTestImag)),'*')
hold on
plot(allLabels - 91, exp(softmaxThetaImag * measurementTestImag))

figure()
plot(signalAngleTest, max(exp(softmaxThetaFull * measurementTestFull)),'*')
hold on
plot(allLabels - 91, exp(softmaxThetaFull * measurementTestFull))
