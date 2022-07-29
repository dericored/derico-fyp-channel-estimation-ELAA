figure('color',[1,1,1]); 
ha=gca;
plot(SNR_dB,nmse_hsamp,'>-','color', '#5F9EA0','linewidth',1.5);
hold on
plot(SNR_dB,nmse_homp,'<-','color', '#A2142F','linewidth',1.5);
hold on
plot(SNR_dB,nmse_LS,'o-','linewidth',1.5, 'color', '#EDB120');
hold on
plot(SNR_dB,nmse_MMSE,'k--','linewidth',1.5);
hold on
grid on
legend('Proposed Hybrid-field SAMP','Hybrid-field OMP','Least Squares','MMSE')
xlabel('SNR (dB)')
ylabel('NMSE (dB)')
ylim([-22 4])
hold off