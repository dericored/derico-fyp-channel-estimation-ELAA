figure('color',[1,1,1]);
ha=gca;
subplot(1,2,1);
plot(step_size_sample,nmse_hsamp,'-o','color', '#5F9EA0','linewidth',1.5);
hold on
grid on
xlabel('Step size {s}')
ylabel('NMSE (dB)')
xlim([1 32])
legend('NMSE of Hybrid-field SAMP')
hold off

subplot(1,2,2);
plot(step_size_sample,exec_time_avg,'-*','color', '#A2142F','linewidth',1.5);
hold on
grid on
xlabel('Step size {s}')
ylabel('Execution time (s)')
xlim([1 32])
legend('Execution time of Hybrid-field SAMP')
hold off