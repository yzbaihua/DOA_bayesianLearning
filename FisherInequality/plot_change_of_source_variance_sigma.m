% this script is used to plot the change of the variance of estiamted
% sigma square

load('MSEResultsSigmaToSource.mat')
load('lowerBoundAlphaResultsSigmaToSource.mat')
load('lowerBoundResultsSigmaToSource.mat')
snrSize = [5,10,15,20,25,30];
data_12_element_10_snapshot_lower_bound = lowerBoundResultsSigmaToSource(:,1);
data_12_element_10_snapshot_lower_bound_alpha = lowerBoundAlphaResultsSigmaToSource(:,1);
data_12_element_10_snapshot_MSE_normal = MSEResultsSigmaToSource(:,1);
data_12_element_10_snapshot_MSE_cluster = MSEResultsSigmaToSource(:,2);

hold on
element12snapshot10Bound = plot(snrSize,data_12_element_10_snapshot_lower_bound,'-*','Color','red');
element12snapshot10BoundAlpha =  plot(snrSize,data_12_element_10_snapshot_lower_bound_alpha,'-o','Color','red');
element12snapshot10MSENormal = plot(snrSize,data_12_element_10_snapshot_MSE_normal,'-*','Color','blue');
element12snapshot10MSECluster =  plot(snrSize,data_12_element_10_snapshot_MSE_cluster,'-o','Color','blue');

legend([element12snapshot10Bound,element12snapshot10BoundAlpha,element12snapshot10MSENormal,element12snapshot10MSECluster],...
    'Lower Bound (13), M = 12, L = 10','Lower Bound with Realization of \alpha (6), M = 12, L = 10', 'MSE, normal array, M = 12, L = 10', 'MSE, cluster array, M = 12, L = 10')
xlabel('SNR')
ylabel('MSE in s')
set(gca, 'YScale', 'log', 'FontSize',12)
grid on