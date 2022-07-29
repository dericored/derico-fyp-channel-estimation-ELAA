clc; 
clear all
close all 

SNR_dB = [0:2:20];
SNR_linear=10.^(SNR_dB./10);
len = length(SNR_linear);
sample = 1000;
sparsity=12;

%%% system parameters
N = 512; % number of beams (transmit antennas)
L = 6; % number of all paths
gamma=0.5; 
Lf = L*gamma; % number of paths for far-field 
Ln = L*(1-gamma); % number of paths for near-field
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
error_LS=zeros(sample,len);
error_MMSE=zeros(sample,len);
energy=zeros(sample,1);

Rh=zeros(N,N);

for s=1:10000
    [h,hf,hn]=generate_hybrid_field_channel(N, Lf, Ln, d, fc,Rmin, Rmax);
    Rh=Rh+h*h';
end
Rh=Rh./(10000);

parfor s=1:sample
    s
    [h,hf,hn] = generate_hybrid_field_channel(N, Lf, Ln, d, fc,Rmin, Rmax);
    
    for iS=1:len
        sigma2=1/SNR_linear(iS);
        P=((rand(M,N)>0.5)*2-1)/sqrt(M); % pilot matrix
        X=eye(N);
        noise = sqrt(sigma2)*(randn(M,1)+1i*randn(M,1))/sqrt(2);
        noise_2 = sqrt(sigma2)*(randn(N,1)+1i*randn(N,1))/sqrt(2);
        y=P*h+noise;
        y_2=X*h+noise_2;
       
        %% the hybrid-field OMP based scheme
        hhat_homp=Hybrid_OMP(y,P,Uf,Un,sparsity*Lf,sparsity*Ln);
        error_homp(s,iS)=sum(abs(hhat_homp-h).^2);

        %% the proposed hybrid-field SAMP based scheme
        step_size=4;
        hhat_hsamp=Hybrid_SAMP(y,P,Uf,Un,step_size,SNR_linear(iS),gamma);
        error_hsamp(s,iS)=sum(abs(hhat_hsamp-h).^2);
        
       %% the LS
       hhat_LS=pinv(X'*X)*X'*y_2;
       error_LS(s,iS)=sum(abs(hhat_LS-h).^2);
       
       %% the MMSE
       hhat_MMSE=Rh*inv(X'*X*Rh+sigma2*eye(N))*X'*y_2;
       error_MMSE(s,iS)=sum(abs(hhat_MMSE-h).^2);
    end
    energy(s)=sum(abs(h).^2);
end
 
nmse_homp = mean(error_homp)/mean(energy);
nmse_hsamp = mean(error_hsamp)/mean(energy);
nmse_LS = mean(error_LS)/mean(energy);
nmse_MMSE = mean(error_MMSE)/mean(energy);

nmse_homp = 10*log10(nmse_homp)
nmse_hsamp = 10*log10(nmse_hsamp)
nmse_LS = 10*log10(nmse_LS)
nmse_MMSE = 10*log10(nmse_MMSE)

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