figure('color',[1,1,1]);
ha=gca;
plot(sparsity_sample,nmse_homp_all(1,:),'<-','linewidth',1.5, 'color', '#EDB120');
hold on
plot(sparsity_sample,nmse_homp_all(2,:),'b>-','linewidth',1.5, 'color', '#5F9EA0');
hold on
plot(sparsity_sample,nmse_homp_all(3,:),'rs-','linewidth',1.5, 'color', '#A2142F');
hold on
plot(sparsity_sample,nmse_homp_all(4,:),'o-','linewidth',1.5, 'color', '#3EB66A');
hold on
grid on
legend('HF-OMP L=6','HF-OMP L=12','HF-OMP L=18','HF-OMP L=24')
xlabel('Assumed Sparsity')
ylabel('NMSE (dB)')
ylim([-12 4])
xlim([1 50])
hold off