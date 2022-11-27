clear variables
SNR =[2,5,10,15,20,25,30];
lowerBoundsLaplacePrior = [0.04036,0.020962,0.006160,0.002007,0.000631,0.000209,0.000062];
lowerBoundsGaussianPrior = [0.019041,0.013006,0.005530,0.0018,0.0006,0.00018,0.00006];
mseLapalce = [1.872386, 1.082846, 0.192330, 0.068541,0.013497,0.001881,0.000626];
mseGaussian = [3.986998, 1.337992, 0.283016, 0.073889,0.018497,0.002196,0.000826];

pLBLaplace = plot(SNR, 10*log10(lowerBoundsLaplacePrior/2),'linewidth',1,'Marker','+','MarkerSize',8);
hold on
pLBGaussian = plot(SNR, 10*log10(lowerBoundsGaussianPrior/2),'linewidth',1,'Marker','^','MarkerSize',8);
mseLaplace = plot(SNR, 10*log10(mseLapalce),'linewidth',1,'Marker','+','MarkerSize',8);
mseGaussian = plot(SNR, 10*log10(mseGaussian),'linewidth',1,'Marker','^','MarkerSize',8);
