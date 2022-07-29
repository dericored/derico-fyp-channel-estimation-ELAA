clc; 
clear all
close all 

sparsity_sample = [1,10:10:50];
len_sparsity = length(sparsity_sample);
sample = 50;
L_sample = [6,12,18,24];
len_L = length(L_sample);
nmse_homp_all = zeros(len_L, len_sparsity);

for l=1:len_L
    %%% system parameters
    N = 512; % number of beams (transmit antennas)
    L = L_sample(l); % number of all paths
    gamma=0.5; 
    Lf = L*gamma; % number of paths for far-field 
    Ln = L*(1-gamma); % number of paths for near-field
    M = 256; % number of pilot overhead

    fc = 30e9; % carrier frequency
    c = 3e8;
    lambda_c = c/fc; % wavelength
    d = lambda_c / 2; % antenna space
    
    SNR_dB = 10;
    SNR_linear = 10^(SNR_dB/10);

    % the far-field angle-domain DFT matrix
    Uf = (1/sqrt(N))*exp(-1i*pi*[0:N-1]'*[-(N-1)/2:1:(N/2)]*(2/N));

    % the near-field polar-domain transform matrix [5]
    Rmin=10;
    Rmax=80;
    eta = 2.5; 
    [Un, label, dict_cell, label_cell] = QuaCode(N, d, lambda_c, eta, Rmin, Rmax);

    error_homp = zeros(sample, len_sparsity);
    energy = zeros(sample, 1);

    parfor s=1:sample
        s
        [h,hf,hn] = generate_hybrid_field_channel(N, Lf, Ln, d, fc,Rmin, Rmax);
        sigma2 = 1/SNR_linear;
        noise = sqrt(sigma2)*(randn(M,1)+1i*randn(M,1))/sqrt(2);
        P=((rand(M,N)>0.5)*2-1)/sqrt(M);
        y=P*h+noise;

        for iS=1:len_sparsity
            % HF-OMP
            sparsity = sparsity_sample(iS);
            hhat_homp = Hybrid_OMP(y,P,Uf,Un,sparsity*Lf,sparsity*Ln);
            error_homp(s,iS)=sum(abs(hhat_homp-h).^2);
        end
        energy(s)=sum(abs(h).^2);
    end

    nmse_homp = mean(error_homp)/mean(energy);
    nmse_homp = 10*log10(nmse_homp);
    
    nmse_homp_all(l,:) = nmse_homp;
end

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