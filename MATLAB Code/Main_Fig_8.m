clc; 
clear all
close all 

SNR_dB = 10;
SNR_linear=10^(SNR_dB/10);
sigma2=1/SNR_linear;
sample = 1000;
sparsity=12;

%%% system parameters
N = 512; % number of beams (transmit antennas)
L_sample = [10:10:50]; % number of all paths
len = length(L_sample);
gamma=0.5;
M = 256; % number of pilot overhead

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

error_homp=zeros(sample,len);
error_hsamp=zeros(sample,len);
energy=zeros(sample,1);

for s=1:sample
    s
    
    for iS=1:len
        L = L_sample(iS);
        Lf = L*gamma;
        Ln = L*(1-gamma);
        [h,hf,hn] = generate_hybrid_field_channel(N, Lf, Ln, d, fc, Rmin, Rmax);
        
        P=((rand(M,N)>0.5)*2-1)/sqrt(M); % pilot matrix
        noise = sqrt(sigma2)*(randn(M,1)+1i*randn(M,1))/sqrt(2);
        y=P*h+noise;
       
        %% the hybrid-field OMP based scheme
        hhat_homp=Hybrid_OMP(y,P,Uf,Un,sparsity*Lf,sparsity*Ln);
        error_homp(s,iS)=sum(abs(hhat_homp-h).^2);

        %% the proposed hybrid-field SAMP based scheme
        step_size=4;
        hhat_hsamp=Hybrid_SAMP(y,P,Uf,Un,step_size,SNR_linear,gamma);
        error_hsamp(s,iS)=sum(abs(hhat_hsamp-h).^2);
        
    end
    energy(s)=sum(abs(h).^2);
end
 
nmse_homp = mean(error_homp)/mean(energy);
nmse_hsamp = mean(error_hsamp)/mean(energy);

nmse_homp = 10*log10(nmse_homp)
nmse_hsamp = 10*log10(nmse_hsamp)

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