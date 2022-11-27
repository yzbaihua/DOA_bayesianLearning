clear variables;
frequency=300e6;
c=3e8;
lambda=c/frequency;
spacing=0.5*lambda;
% signal to noise ratio, we define it by SNR = max(signal)/noise
SNR=2; 

% signal=[1.1*exp(1i*1.05);0.8*exp(1i*0.3);1.5*exp(1i*1.3)]; % the signal 
signal=[sqrt(2)*exp(1i*pi/4);sqrt(2)*exp(1i*pi/4)]; % the signal 

% signalVector=[real(signal);imag(signal)];
deltaDis = 103*lambda;
loc =[1,2,3, [1,2,3]+deltaDis, [1,2,3]+2*deltaDis,[1,2,3]+3*deltaDis];

signalAngleSet = combnk(-90:90,2);
trainSize = length(signalAngleSet);
trainMatrixRealPart = zeros(length(loc),trainSize);
trainMatrixImagPart = zeros(length(loc),trainSize);
trainMatrixFull = zeros(2*length(loc),trainSize);
trainLabel = zeros(181,trainSize);
for trainIdx = 1 : trainSize
%     signalAngle=randi([-90,90],1,1);
    signalAngle= signalAngleSet(trainIdx,:);
%     signal=sqrt(2*rand)*exp(1i*pi/(4*rand)); % the signal
    signal=(sqrt(2)*exp(1i*pi/(4)))*ones(length(signalAngle),1); % the signal
    trainLabel(signalAngle + 91,trainIdx) = 1;
    actualReceivedMatrix = exp(-1i*loc'*2*pi*sind(signalAngle)/lambda)*signal;
    trainMatrixRealPart(:,trainIdx) = real(actualReceivedMatrix);
    trainMatrixImagPart(:,trainIdx) = imag(actualReceivedMatrix); 
    trainMatrixFull(:,trainIdx) = [real(actualReceivedMatrix);imag(actualReceivedMatrix)];
end

inputSize = length(loc);
numCls = size(trainLabel,1);
softmaxLambda = 1e-4;
options.maxIter = 500;

signalAngleTest = [50,-20];
signalTest = [sqrt(2)*exp(1i*pi/4);sqrt(2)*exp(1i*pi/4)];%sqrt(2)*exp(1i*pi/4)];
noiseTest = my_awgn(signalTest,SNR)-signalTest;
% noiseless
measurementTest = exp(-1i*loc'*2*pi*sind(signalAngleTest)/lambda)...
    *(signalTest);
% % noisy
% measurementTest = exp(-1i*loc'*2*pi*sind(signalAngleTest)/lambda)...
%     *(signalTest+noiseTest);

% training part ---------
% predict real part
softmaxModelReal = softmaxTrainMultipleSource(inputSize,numCls,softmaxLambda,...
    trainMatrixRealPart,trainLabel,options);
saeSoftmaxOptThetaReal = softmaxModelReal.optTheta(:);
softmaxThetaReal = reshape(saeSoftmaxOptThetaReal(1:inputSize*numCls),numCls,inputSize);

% predict imag part
softmaxModelImag = softmaxTrainMultipleSource(inputSize,numCls,softmaxLambda,...
    trainMatrixImagPart,trainLabel,options);
saeSoftmaxOptThetaImag = softmaxModelImag.optTheta(:);
softmaxThetaImag = reshape(saeSoftmaxOptThetaImag(1:inputSize*numCls),numCls,inputSize);

% predict full part
softmaxModelFull = softmaxTrainMultipleSource(2*inputSize,numCls,softmaxLambda,...
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

% accuracy test
% result = zeros(181,1);
% for i  = 1 : 181
%     signalAngleTest = i - 91;
%     signalTest = [sqrt(2)*exp(1i*pi/4)];%sqrt(2)*exp(1i*pi/4);sqrt(2)*exp(1i*pi/4)];
%     noiseTest = my_awgn(signalTest,SNR)-signalTest;
%     % noiseless
%     measurementTestNoiseless = exp(-1i*loc'*2*pi*sind(signalAngleTest)/lambda)...
%         *(signalTest);
% 
%     % noisy_1 or 2
%     % % noisy_1
%     measurementTest = exp(-1i*loc'*2*pi*sind(signalAngleTest)/lambda)...
%         *(signalTest+noiseTest);
% 
%     measurementTestReal = real(measurementTest);
%     measurementTestImag = imag(measurementTest);
%     measurementTestFull = [measurementTestReal;measurementTestImag];
%     probability = exp(softmaxThetaFull * measurementTestFull)./sum(exp(softmaxThetaFull * measurementTestFull));
%     [~,result(i)] = max(probability);
% end
% end of accuracy test 

figure()
plot(signalAngleTest, ones(length(signalAngleTest),1),'*')
hold on
plot(-90:90, exp(softmaxThetaReal * measurementTestReal))
figure()
plot(signalAngleTest, max(exp(softmaxThetaImag * measurementTestImag)),'*')
hold on
plot(-90:90, exp(softmaxThetaImag * measurementTestImag))

figure()
probability = exp(softmaxThetaFull * measurementTestFull)./sum(exp(softmaxThetaFull * measurementTestFull));
h1 = plot(signalAngleTest, max(probability),'*','MarkerSize',10,'color','red');
hold on
h2 = plot(-90:90,probability,'color','blue');
ylabel('Probability(%)')
xlabel('Theta ($\theta$)','Interpreter','latex')
lg = legend(h1,'Actual direction');
set(gca,'fontsize',16)
