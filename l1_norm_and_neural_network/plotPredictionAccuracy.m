clear variable
% the accuracy under 2db, 5db, 10db, 15dB 20db 25 db and 30 db
oneSoueceAccuracy = [42.3 62.3 75 82.3 86.6 87.8 89.1];
oneSoueceAccuracyNeuralwork = [56.8 66.5 79.4 91.7 94.1 95.8 97.3];

twoSoueceAccuracy = [32.8 35.75 39.8 43.1 45.15 45.8 46.65]; 
twoSoueceAccuracyNeuralwork = [26.75 36.5 40.15 44.2 47.85 48 48.5];

SNR = [2 5 10 15 20 25 30];
plot(SNR, oneSoueceAccuracy, '-^', 'MarkerSize', 10)
xlabel('SNR(dB)', 'FontSize', 10)
ylabel('Prediction Accuracy(%)', 'FontSize', 10)
grid on

figure()
plot(SNR, twoSoueceAccuracy, '-^', 'MarkerSize', 10)
xlabel('SNR(dB)', 'FontSize', 10)
ylabel('Prediction Accuracy(%)', 'FontSize', 10)
grid on

figure()
onePlot = plot(SNR, oneSoueceAccuracy, '-^', 'MarkerSize', 10);
hold on
twoPlot = plot(SNR, oneSoueceAccuracyNeuralwork, '-^', 'MarkerSize', 10, 'Color', 'red');
legend([onePlot twoPlot],{'Softmax Classifier','Single-layer Network'},'FontSize', 12)
grid on
xlabel('SNR(dB)', 'FontSize', 10)
ylabel('Prediction Accuracy(%)', 'FontSize', 10)

figure()
onePlot = plot(SNR, twoSoueceAccuracy, '-^', 'MarkerSize', 10);
hold on
twoPlot = plot(SNR, twoSoueceAccuracyNeuralwork, '-^', 'MarkerSize', 10, 'Color', 'red');
legend([onePlot twoPlot],{'Softmax Classifier','Single-layer Network'},'FontSize', 12)
grid on
xlabel('SNR(dB)', 'FontSize', 10)
ylabel('Prediction Accuracy(%)', 'FontSize', 10)