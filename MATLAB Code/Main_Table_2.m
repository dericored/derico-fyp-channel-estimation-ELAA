clc; 
clear all
close all 

SNR_dB = 5;
SNR_linear=10^(SNR_dB/10);
sigma2=1/SNR_linear;
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

error_homp=zeros(sample,1);
error_hsamp=zeros(sample,1);
error_LS=zeros(sample,1);
error_MMSE=zeros(sample,1);
energy=zeros(sample,1);

time_homp=zeros(sample,1);
time_hsamp=zeros(sample,1);
time_LS=zeros(sample,1);
time_MMSE=zeros(sample,1);

Rh=zeros(N,N);

for s=1:10000
    [h,hf,hn]=generate_hybrid_field_channel(N, Lf, Ln, d, fc,Rmin, Rmax);
    Rh=Rh+h*h';
end
Rh=Rh./(10000);

parfor s=1:sample
    s
    [h,hf,hn] = generate_hybrid_field_channel(N, Lf, Ln, d, fc,Rmin, Rmax);
    
    P=((rand(M,N)>0.5)*2-1)/sqrt(M); % pilot matrix
    X=eye(N);
    noise = sqrt(sigma2)*(randn(M,1)+1i*randn(M,1))/sqrt(2);
    noise_2 = sqrt(sigma2)*(randn(N,1)+1i*randn(N,1))/sqrt(2);
    y=P*h+noise;
    y_2=X*h+noise_2;

    %% the hybrid-field OMP based scheme
    tic;
    hhat_homp=Hybrid_OMP(y,P,Uf,Un,sparsity*Lf,sparsity*Ln);
    t=toc;
    error_homp(s)=sum(abs(hhat_homp-h).^2);
    time_homp(s)=t;

    %% the proposed hybrid-field SAMP based scheme
    step_size=4;
    tic;
    hhat_hsamp=Hybrid_SAMP(y,P,Uf,Un,step_size,SNR_linear,gamma);
    t=toc;
    error_hsamp(s)=sum(abs(hhat_hsamp-h).^2);
    time_hsamp(s)=t;

   %% the LS
   tic;
   hhat_LS=pinv(X'*X)*X'*y_2;
   t=toc;
   error_LS(s)=sum(abs(hhat_LS-h).^2);
   time_LS(s)=t;

   %% the MMSE
   tic;
   hhat_MMSE=Rh*inv(X'*X*Rh+sigma2*eye(N))*X'*y_2;
   t=toc;
   error_MMSE(s)=sum(abs(hhat_MMSE-h).^2);
   time_MMSE(s)=t;
   
   energy(s)=sum(abs(h).^2);
end
 
nmse_homp = mean(error_homp)/mean(energy);
nmse_hsamp = mean(error_hsamp)/mean(energy);
nmse_LS = mean(error_LS)/mean(energy);
nmse_MMSE = mean(error_MMSE)/mean(energy);

nmse_homp = 10*log10(nmse_homp);
nmse_hsamp = 10*log10(nmse_hsamp);
nmse_LS = 10*log10(nmse_LS);
nmse_MMSE = 10*log10(nmse_MMSE);

time_homp_avg = mean(time_homp);
time_hsamp_avg = mean(time_hsamp);
time_LS_avg = mean(time_LS);
time_MMSE_avg = mean(time_MMSE);

label = ["HF-OMP"; "HF-SAMP"; "Least Squares"; "MMSE"];
nmse = [nmse_homp; nmse_hsamp; nmse_LS; nmse_MMSE];
time = [time_homp_avg; time_hsamp_avg; time_LS_avg; time_MMSE_avg];

performance = table(label, nmse, time)