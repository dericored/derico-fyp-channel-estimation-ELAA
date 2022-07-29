clc; 
clear all
close all 

sparsity_sample = [1,10:5:50];
len_sparsity = length(sparsity_sample);
sample = 1000;
SNR_dB = [5,10];
SNR_linear = 10.^(SNR_dB./10);
len_SNR = length(SNR_dB);

%%% system parameters
N = 512; % number of beams (transmit antennas)
L = 6;
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

error_homp = zeros(sample, len_sparsity);
error_hsamp = zeros(sample, len_sparsity);
error_LS = zeros(sample, len_sparsity);
error_MMSE = zeros(sample, len_sparsity);
energy = zeros(sample, 1);

Rh=zeros(N,N);

for s=1:10000
    [h,hf,hn]=generate_hybrid_field_channel(N, Lf, Ln, d, fc,Rmin, Rmax);
    Rh=Rh+h*h';
end
Rh=Rh./(10000);

nmse_homp = {};
nmse_hsamp = {};
nmse_LS = {};
nmse_MMSE = {};

for n=1:len_SNR
    SNR = SNR_linear(n);
    for s=1:sample
        s
        [h,hf,hn] = generate_hybrid_field_channel(N, Lf, Ln, d, fc,Rmin, Rmax);
        sigma2 = 1/SNR;
        noise = sqrt(sigma2)*(randn(M,1)+1i*randn(M,1))/sqrt(2);
        noise_2 = sqrt(sigma2)*(randn(N,1)+1i*randn(N,1))/sqrt(2);
        P=((rand(M,N)>0.5)*2-1)/sqrt(M);
        X=eye(N,N);
        y=P*h+noise;
        y_2=X*h+noise_2;

        for iS=1:len_sparsity
           %% HF-OMP
           sparsity = sparsity_sample(iS);
           hhat_homp = Hybrid_OMP(y,P,Uf,Un,sparsity*Lf,sparsity*Ln);
           error_homp(s,iS)=sum(abs(hhat_homp-h).^2);

           %% HF-SAMP
           step_size=4;
           hhat_hsamp=Hybrid_SAMP(y,P,Uf,Un,step_size,SNR,gamma);
           error_hsamp(s,iS)=sum(abs(hhat_hsamp-h).^2);

           %% the LS
           hhat_LS=pinv(X'*X)*X'*y_2;
           error_LS(s,iS)=sum(abs(hhat_LS-h).^2);;

           %% the MMSE
           hhat_MMSE=Rh*inv(X'*X*Rh+sigma2*eye(N))*X'*y_2;
           error_MMSE(s,iS)=sum(abs(hhat_MMSE-h).^2);
        end
        energy(s)=sum(abs(h).^2);
    end

    nmse_homp_temp = mean(error_homp)/mean(energy);
    nmse_hsamp_temp = mean(error_hsamp)/mean(energy);
    nmse_LS_temp = mean(error_LS)/mean(energy);
    nmse_MMSE_temp = mean(error_MMSE)/mean(energy);

    nmse_homp{n} = 10*log10(nmse_homp_temp);
    nmse_hsamp{n} = 10*log10(nmse_hsamp_temp);
    nmse_LS{n} = 10*log10(nmse_LS_temp);
    nmse_MMSE{n} = 10*log10(nmse_MMSE_temp);

end

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

