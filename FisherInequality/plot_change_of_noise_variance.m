% this script is used to plot the change of the variance of estiamted
% sigma square

load('lowerBoundResults.mat')
load('MSEresults.mat')
snrSize = [5,10,15,20,25,30];
data_6_element_10_snapshot_lower_bound = lowerBoundResults(:,1);
data_6_element_10_snapshot_MSE = MSEresults(:,1);
data_9_element_10_snapshot_lower_bound = lowerBoundResults(:,2);
data_9_element_10_snapshot_MSE = MSEresults(:,2);
data_12_element_10_snapshot_lower_bound = lowerBoundResults(:,3);
data_12_element_10_snapshot_MSE = MSEresults(:,3);
data_15_element_10_snapshot_lower_bound = lowerBoundResults(:,4);
data_15_element_10_snapshot_MSE = MSEresults(:,4);
data_18_element_10_snapshot_lower_bound = lowerBoundResults(:,5);
data_18_element_10_snapshot_MSE = MSEresults(:,5);

hold on
element6snapshot10Bound = plot(snrSize,10*log10(data_6_element_10_snapshot_lower_bound),'-*','Color','red');
element6snapshot10MSE =  plot(snrSize,10*log10(data_6_element_10_snapshot_MSE),'-*','Color','blue');
% element9snapshot10Bound = plot(snrSize,data_9_element_10_snapshot_lower_bound,'-o','Color','red');
% element9snapshot10MSE =  plot(snrSize,data_9_element_10_snapshot_MSE,'-o','Color','blue');
element12snapshot10Bound = plot(snrSize,10*log10(data_12_element_10_snapshot_lower_bound),'-+','Color','red');
element12snapshot10MSE =  plot(snrSize,10*log10(data_12_element_10_snapshot_MSE),'-+','Color','blue');
% element15snapshot10Bound = plot(snrSize,data_15_element_10_snapshot_lower_bound,'-x','Color','red');
% element15snapshot10MSE =  plot(snrSize,data_15_element_10_snapshot_MSE,'-x','Color','blue');
element18snapshot10Bound = plot(snrSize,10*log10(data_18_element_10_snapshot_lower_bound),'-d','Color','red');
element18snapshot10MSE =  plot(snrSize,10*log10(data_18_element_10_snapshot_MSE),'-d','Color','blue');
legend([element6snapshot10Bound,element12snapshot10Bound,element18snapshot10Bound,element6snapshot10MSE,element12snapshot10MSE,element18snapshot10MSE],...
    'Lower Bound, M = 6','Lower Bound, M = 12', 'Lower Bound, M = 18', 'MSE, M = 6','MSE, M = 12', 'MSE, M = 18')
xlabel('SNR')
ylabel('MSE in \sigma^{2}')
set(gca, 'YScale', 'log', 'FontSize',12)
grid on

figure()
snrSize = [5,10,15,20,25,30];
data_6_element_50_snapshot_lower_bound = lowerBoundResults(:,6);
data_6_element_50_snapshot_MSE = MSEresults(:,6);
data_9_element_50_snapshot_lower_bound = lowerBoundResults(:,7);
data_9_element_50_snapshot_MSE = MSEresults(:,7);
data_12_element_50_snapshot_lower_bound = lowerBoundResults(:,8);
data_12_element_50_snapshot_MSE = MSEresults(:,8);
data_15_element_50_snapshot_lower_bound = lowerBoundResults(:,9);
data_15_element_50_snapshot_MSE = MSEresults(:,9);
data_18_element_50_snapshot_lower_bound = lowerBoundResults(:,10);
data_18_element_50_snapshot_MSE = MSEresults(:,10);

hold on
element6snapshot50Bound = plot(snrSize,data_6_element_50_snapshot_lower_bound,'-*','Color','red');
element6snapshot50MSE =  plot(snrSize,data_6_element_50_snapshot_MSE,'-*','Color','blue');
% element9snapshot50Bound = plot(snrSize,data_9_element_50_snapshot_lower_bound,'-o','Color','red');
% element9snapshot50MSE =  plot(snrSize,data_9_element_50_snapshot_MSE,'-o','Color','blue');
element12snapshot50Bound = plot(snrSize,data_12_element_50_snapshot_lower_bound,'-+','Color','red');
element12snapshot50MSE =  plot(snrSize,data_12_element_50_snapshot_MSE,'-+','Color','blue');
% element15snapshot50Bound = plot(snrSize,data_15_element_50_snapshot_lower_bound,'-x','Color','red');
% element15snapshot50MSE =  plot(snrSize,data_15_element_50_snapshot_MSE,'-x','Color','blue');
element18snapshot50Bound = plot(snrSize,data_18_element_50_snapshot_lower_bound,'-d','Color','red');
element18snapshot50MSE =  plot(snrSize,data_18_element_50_snapshot_MSE,'-d','Color','blue');

legend([element6snapshot50Bound,element12snapshot50Bound,element18snapshot50Bound,element6snapshot50MSE,element12snapshot50MSE,element18snapshot50MSE],...
    'Lower Bound, M = 6','Lower Bound, M = 12', 'Lower Bound, M = 18', 'MSE, M = 6','MSE, M = 12', 'MSE, M = 18')
xlabel('SNR')
ylabel('MSE in \sigma^{2}')
set(gca, 'YScale', 'log', 'FontSize',12)
grid on



