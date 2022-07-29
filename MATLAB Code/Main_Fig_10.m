clc; 
clear all
close all 

%%% system parameters
N = 512; % number of beams (transmit antennas)
L = 12; % number of all paths
gamma=0.5; 
Lf = L*gamma; % number of paths for far-field 
Ln = L*(1-gamma); % number of paths for near-field
M = 256; % number of pilot overhead

SNR_dB = 10;
SNR_linear=10^(SNR_dB/10);
sigma2 = 1/SNR_linear;
step_size_sample = [1:1:32];
len = length(step_size_sample);
sample = 1000;

fc = 30e9; % carrier frequency
c = 3e8;
lambda_c = c/fc; % wavelength 
d = lambda_c / 2; % antenna space

% the far-field angle-domain DFT matrix
Uf = (1/sqrt(N))*exp(-1i*pi*[0:N-1]'*[-(N-1)/2:1:(N/2)]*(2/N));

% the near-field polar-domain transform matrix [5]
Rmin=10;
Rmax=80;
eta = 2.5; 
[Un, label, dict_cell, label_cell] = QuaCode(N, d, lambda_c, eta, Rmin, Rmax);

error_hsamp=zeros(sample,len);
exec_time=zeros(sample,len);

parfor s=1:sample
    s
    [h,hf,hn] = generate_hybrid_field_channel(N, Lf, Ln, d, fc,Rmin, Rmax);
    
    for iS=1:len
        P=((rand(M,N)>0.5)*2-1)/sqrt(M); % pilot matrix
        noise = sqrt(sigma2)*(randn(M,1)+1i*randn(M,1))/sqrt(2);
        y=P*h+noise;
       
        %% the proposed hybrid-field SAMP based scheme
        step_size=step_size_sample(iS);
        tic;
        hhat_hsamp=Hybrid_SAMP(y,P,Uf,Un,step_size,SNR_linear,gamma);
        t=toc;
        error_hsamp(s,iS)=sum(abs(hhat_hsamp-h).^2);
        exec_time(s,iS)=t;
    end
    energy(s)=sum(abs(h).^2);
end
 
nmse_hsamp = mean(error_hsamp)/mean(energy);
nmse_hsamp = 10*log10(nmse_hsamp)

exec_time_avg = mean(exec_time);

performance=table(step_size_sample', nmse_hsamp', exec_time_avg')

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