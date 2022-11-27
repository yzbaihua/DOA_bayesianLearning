% plot logistic function
z = -5:5;
g_z = 1./(1+exp(-1*z));
plot(z,g_z,'linewidth',2);
grid on
xlabel('z','fontsize',15);
ylabel('g(z)','fontsize',15);
% set(gca,'fontsize',15)