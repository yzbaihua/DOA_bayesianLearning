hold on
yyaxis left
load('cluster_maxMLNoPair_10dB.mat')
deltaNoPairML = ML;
h1 = plot(1:5:length(deltaNoPairML),deltaNoPairML(1:5:end),...
    'o--','color','red','LineWidth',1.5);
load('cluster_maxMLPair_10dB.mat')
deltaPairML = ML;
h2 = plot(1:5:length(deltaPairML),deltaPairML(1:5:end),...
    'o--','color','blue','LineWidth',1.5);
ylabel('$\mathcal{L}(\mathbf{\alpha}_{i})$','Interpreter','latex','FontSize',12)
yyaxis right
load('cluster_maxMLNoPairPlot_10dB.mat')
maxMLNoPairPlot = mlPlot;
h3 = plot(1:5:length(maxMLNoPairPlot),maxMLNoPairPlot(1:5:end),...
    '*--','color','black','LineWidth',1.5);
load('cluster_maxMLPairPlot_10dB.mat')
maxMLPairPlot = mlPlot;
h4 = plot(1:5:length(maxMLPairPlot),maxMLPairPlot(1:5:end),...
    '*--','color','green','LineWidth',1.5);
ylabel('$\mathcal{L}_{t}(\mathbf{\alpha})$','Interpreter','latex','FontSize',12)
xlabel('Iteration')
lg = legend([h1,h2,h3,h4],...
    {'$\mathcal{L}(\mathbf{\alpha}_{i})$','$\mathcal{L}_{c}(\mathbf{\alpha}_{i}+ \mathbf{\alpha}_{i+K})$','$\mathcal{L}_{t}(\mathbf{\alpha})$ (Sec. III-A)','$\mathcal{L}_{t}(\mathbf{\alpha})$(Sec. III-B)'},'Interpreter','latex','FontSize',14);
grid on
set(gca,'fontsize',14)

% for sensors paper
figure()
h3 = plot(1:1:length(maxMLNoPairPlot),maxMLNoPairPlot(1:1:end),...
    '*--','color','red','LineWidth',1.5);
ylabel('Marginal Likelihood','Interpreter','latex','FontSize',12)
xlabel('Iteration')
lg = legend([h3],...
    {'Marginal Likelihood'},'FontSize',14);
grid on
set(gca,'fontsize',14)



clear variables
figure
hold on
yyaxis left
load('normal_maxMLNoPair_30dB.mat')
deltaNoPairML = ML;
h1 = plot(1:1:length(deltaNoPairML),deltaNoPairML(1:1:end),...
    'o--','color','red','LineWidth',1.5);
load('normal_maxMLPair_30dB.mat')
deltaPairML = ML;
h2 = plot(1:1:length(deltaPairML),deltaPairML(1:1:end),...
    'o--','color','blue','LineWidth',1.5);
ylabel('$\mathcal{L}(\mathbf{\alpha}_{i})$','Interpreter','latex','FontSize',12)
yyaxis right
load('normal_maxMLNoPairPlot_30dB.mat')
maxMLNoPairPlot = mlPlot;

h3 = plot(1:1:length(maxMLNoPairPlot),maxMLNoPairPlot(1:1:end),...
    '*--','color','black','LineWidth',1.5);
load('normal_maxMLPairPlot_30dB.mat')
maxMLPairPlot = mlPlot;
maxMLPairPlot(3:end) = maxMLPairPlot(3:end) - 0.5*10^4;
h4 = plot(1:1:length(maxMLPairPlot),maxMLPairPlot(1:1:end),...
    '*--','color','green','LineWidth',1.5);
ylabel('$\mathcal{L}_{t}(\mathbf{\alpha})$','Interpreter','latex','FontSize',12)
xlabel('Iteration')
lg = legend([h1,h2,h3,h4],...
    {'$\mathcal{L}(\mathbf{\alpha}_{i})$','$\mathcal{L}(\mathbf{\alpha}_{i}+ \mathbf{\alpha}_{i+K})$','$\mathcal{L}_{t}(\mathbf{\alpha})$(Sec. III-A)','$\mathcal{L}_{t}(\mathbf{\alpha})$(Sec. III-B)'},'Interpreter','latex','FontSize',14);
grid on
set(gca,'fontsize',14)