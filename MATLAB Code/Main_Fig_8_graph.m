figure('color',[1,1,1]); 
ha=gca;
plot(L_sample,nmse_homp,'rs-','color', '#A2142F','linewidth',1.5);
hold on
plot(L_sample,nmse_hsamp,'b>-','color', '#5F9EA0','linewidth',1.5);
hold on
grid on
legend('Hybrid-field OMP','Proposed Hybrid-field SAMP')
xlabel('Number of scatter paths L')
ylabel('NMSE (dB)')
ylim([-12 4])
xticks([10 20 30 40 50])
hold off