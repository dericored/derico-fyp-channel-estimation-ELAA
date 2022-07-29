% SNR = 5 dB
figure(1); 
ha=gca;
plot(sparsity_sample,nmse_hsamp{1},'>-','color', '#5F9EA0','linewidth',1.5);
hold on
plot(sparsity_sample,nmse_homp{1},'<-','color', '#A2142F','linewidth',1.5);
hold on
plot(sparsity_sample,nmse_LS{1},'o-','linewidth',1.5, 'color', '#EDB120');
hold on
plot(sparsity_sample,nmse_MMSE{1},'k--','linewidth',1.5);
hold on
grid on
legend('Proposed Hybrid-field SAMP','Hybrid-field OMP','Least Squares','MMSE')
xlabel('Assumed sparsity')
ylabel('NMSE (dB)')
ylim([-10 2])
xlim([1 50])
hold off

% SNR = 10 dB
figure(2); 
ha=gca;
plot(sparsity_sample,nmse_hsamp{2},'>-','color', '#5F9EA0','linewidth',1.5);
hold on
plot(sparsity_sample,nmse_homp{2},'<-','color', '#A2142F','linewidth',1.5);
hold on
plot(sparsity_sample,nmse_LS{2},'o-','linewidth',1.5, 'color', '#EDB120');
hold on
plot(sparsity_sample,nmse_MMSE{2},'k--','linewidth',1.5);
hold on
grid on
legend('Proposed Hybrid-field SAMP','Hybrid-field OMP','Least Squares','MMSE')
xlabel('Assumed sparsity')
ylabel('NMSE (dB)')
ylim([-14 0])
xlim([1 50])
hold off