% This Matlab script can be used to generate Fig. 10 in the paper:
% R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

clear;
clc;

M = 6;  %% number of transmit antennas
N = 100; %% number of reflecting elements
K = 3;  %% number of primary receivers
phi = pi/4;  %% QPSK for PRs
global sigma2
sigma2 = 10^-11; %% noise power -80dBm

dist_SR = 0.5;  %%
dist_PR = 2.5;  %% QoS of PRs
dist = [dist_PR.*ones(1,K) dist_SR];
omega = 4.*ones(1,K);


N_sim = 1000;  %% number of simulations
power_range = (22:1:28);

SER_my_SR_L1 = zeros(1,length(power_range));
SER_my_2_SR_L1 = zeros(1,length(power_range));
SER_my_4_SR_L1 = zeros(1,length(power_range));

SER_my_SR_L2 = zeros(1,length(power_range));
SER_my_2_SR_L2 = zeros(1,length(power_range));
SER_my_4_SR_L2 = zeros(1,length(power_range));

SER_my_SR_L3 = zeros(1,length(power_range));
SER_my_2_SR_L3 = zeros(1,length(power_range));
SER_my_4_SR_L3 = zeros(1,length(power_range));

B2 = 2;
B4 = 4;

d_ar = 10;  %% BS-IRS distance
d_rs = 20; %%%%
d_ru = 100;

belta1 = sqrt(10^(0.3)/(1+10^(0.3)));
belta2 = sqrt(1/(1+10^(0.3)));

H_au = zeros(K+1,M);
H_ru = zeros(K+1,N);

Nmax = 9;
res_th = 1e-3;

L1 = 8;
L2 = 12;
L3 = 16;
Ns = 9600;

R = zeros(K+1,Ns);
R2 = zeros(K+1,Ns);
R4 = zeros(K+1,Ns);
for sim = 1:N_sim
    sim
    tic
    d_as = sqrt(d_ar^2+d_rs^2)+(d_ar+d_rs-sqrt(d_ar^2+d_rs^2))*rand;
    d_au = sqrt(d_ar^2+d_ru^2)+(d_ar+d_ru-sqrt(d_ar^2+d_ru^2)).*rand(1,K);

    H_ar = sqrt(10^(-3)*d_ar^(-2.5))*(belta1*channel_ar(M,N)+belta2*(randn(N,M)+1i*randn(N,M))/sqrt(2));
    H_ru(1:K,:) = sqrt(10^(-3)*d_ru^(-3))*(randn(K,N)+1i*randn(K,N))/sqrt(2);
    H_ru(K+1,:) = sqrt(10^(-3)*d_rs^(-3))*(randn(1,N)+1i*randn(1,N))/sqrt(2);
    for i = 1:1:K
        H_au(i,:) = sqrt(10^(-3)*d_au(i)^(-3))*(randn(1,M)+1i*randn(1,M))/sqrt(2);
    end
    H_au(K+1,:) = sqrt(10^(-3)*d_as^(-3))*(randn(1,M)+1i*randn(1,M))/sqrt(2);
    
    Noise = sqrt(0.5*sigma2)*(randn(K+1,Ns)+1i*randn(K+1,Ns));
    S_PR = exp(1i*2*pi*rand(K,Ns));
    S_SR_L1 = exp(1i*2*pi*rand(1,Ns/L1));
    temp = reshape(repmat(S_SR_L1.',1,L1).',Ns,1);
    [S_BS_L1,index_BS_L1] = get_adaptive_modulate([S_PR;temp.'],[omega 2]);
    S_SR_L2 = S_SR_L1(1:Ns/L2);
    temp = reshape(repmat(S_SR_L2.',1,L2).',Ns,1);
    [S_BS_L2,index_BS_L2] = get_adaptive_modulate([S_PR;temp.'],[omega 2]);
    S_SR_L3 = S_SR_L1(1:Ns/L3);
    temp = reshape(repmat(S_SR_L3.',1,L3).',Ns,1);
    [S_BS_L3,index_BS_L3] = get_adaptive_modulate([S_PR;temp.'],[omega 2]);

    [S_PR,index_PR] = get_adaptive_modulate(S_PR,omega);  %% modulate PRs' information

    [S_SR_L1,index_SR_L1] = get_adaptive_modulate(S_SR_L1,2);  %% modulate SR's information
    index_theta_L1 = reshape(repmat(index_SR_L1,L1,1),Ns,1)';
    index_SR_L2 = index_SR_L1(1:Ns/L2);
    index_theta_L2 = reshape(repmat(index_SR_L2,L2,1),Ns,1)';
    index_SR_L3 = index_SR_L1(1:Ns/L3);
    index_theta_L3 = reshape(repmat(index_SR_L3,L3,1),Ns,1)';
    S_SR_L2 = S_SR_L1(1:Ns/L2);
    S_SR_L3 = S_SR_L1(1:Ns/L3);

    for power_index = 1:length(power_range)
        power = power_range(power_index)

        [X_my,Theta_my,t] = get_X_theta(H_au,H_ar,H_ru,dist,power,Nmax,res_th);
        [X_my2,Theta_my2,t2] = get_X_theta_b(H_au,H_ar,H_ru,dist,power,B2,Nmax,res_th);
        [X_my4,Theta_my4,t4] = get_X_theta_b(H_au,H_ar,H_ru,dist,power,B4,Nmax,res_th);

        for i = 1:1:Ns
            R(:,i) = (H_au + H_ru*diag(Theta_my(:,index_theta_L1(i)))*H_ar)*X_my(:,index_PR(i))+ Noise(:,i);
            R2(:,i) = (H_au + H_ru*diag(Theta_my2(:,index_theta_L1(i)))*H_ar)*X_my2(:,index_PR(i))+ Noise(:,i);
            R4(:,i) = (H_au + H_ru*diag(Theta_my4(:,index_theta_L1(i)))*H_ar)*X_my4(:,index_PR(i))+ Noise(:,i);
        end
        SER_my_SR_L1(power_index) = SER_my_SR_L1(power_index) + get_SER(sum(reshape(R(K+1,:),L1,Ns/L1))/L1,S_SR_L1,2);
        SER_my_2_SR_L1(power_index) = SER_my_2_SR_L1(power_index) + get_SER(sum(reshape(R2(K+1,:),L1,Ns/L1))/L1,S_SR_L1,2);
        SER_my_4_SR_L1(power_index) = SER_my_4_SR_L1(power_index) + get_SER(sum(reshape(R4(K+1,:),L1,Ns/L1))/L1,S_SR_L1,2);

        for i = 1:1:Ns
            R(K+1,i) = (H_au(K+1,:) + H_ru(K+1,:)*diag(Theta_my(:,index_theta_L2(i)))*H_ar)*X_my(:,index_PR(i))+ Noise(K+1,i);
            R2(K+1,i) = (H_au(K+1,:) + H_ru(K+1,:)*diag(Theta_my2(:,index_theta_L2(i)))*H_ar)*X_my2(:,index_PR(i))+ Noise(K+1,i);
            R4(K+1,i) = (H_au(K+1,:) + H_ru(K+1,:)*diag(Theta_my4(:,index_theta_L2(i)))*H_ar)*X_my4(:,index_PR(i))+ Noise(K+1,i);
        end
        SER_my_SR_L2(power_index) = SER_my_SR_L2(power_index) + get_SER(sum(reshape(R(K+1,:),L2,Ns/L2))/L2,S_SR_L2,2);
        SER_my_2_SR_L2(power_index) = SER_my_2_SR_L2(power_index) + get_SER(sum(reshape(R2(K+1,:),L2,Ns/L2))/L2,S_SR_L2,2);
        SER_my_4_SR_L2(power_index) = SER_my_4_SR_L2(power_index) + get_SER(sum(reshape(R4(K+1,:),L2,Ns/L2))/L2,S_SR_L2,2);

        for i = 1:1:Ns
            R(K+1,i) = (H_au(K+1,:) + H_ru(K+1,:)*diag(Theta_my(:,index_theta_L3(i)))*H_ar)*X_my(:,index_PR(i))+ Noise(K+1,i);
            R2(K+1,i) = (H_au(K+1,:) + H_ru(K+1,:)*diag(Theta_my2(:,index_theta_L3(i)))*H_ar)*X_my2(:,index_PR(i))+ Noise(K+1,i);
            R4(K+1,i) = (H_au(K+1,:) + H_ru(K+1,:)*diag(Theta_my4(:,index_theta_L3(i)))*H_ar)*X_my4(:,index_PR(i))+ Noise(K+1,i);
        end
        SER_my_SR_L3(power_index) = SER_my_SR_L3(power_index) + get_SER(sum(reshape(R(K+1,:),L3,Ns/L3))/L3,S_SR_L3,2);
        SER_my_2_SR_L3(power_index) = SER_my_2_SR_L3(power_index) + get_SER(sum(reshape(R2(K+1,:),L3,Ns/L3))/L3,S_SR_L3,2);
        SER_my_4_SR_L3(power_index) = SER_my_4_SR_L3(power_index) + get_SER(sum(reshape(R4(K+1,:),L3,Ns/L3))/L3,S_SR_L3,2);

    end
    toc
end

SER_my_SR_L1 = SER_my_SR_L1/sim;
SER_my_2_SR_L1 = SER_my_2_SR_L1/sim;
SER_my_4_SR_L1 = SER_my_4_SR_L1/sim;

SER_my_SR_L2 = SER_my_SR_L2/sim;
SER_my_2_SR_L2 = SER_my_2_SR_L2/sim;
SER_my_4_SR_L2 = SER_my_4_SR_L2/sim;

SER_my_SR_L3 = SER_my_SR_L3/sim;
SER_my_2_SR_L3 = SER_my_2_SR_L3/sim;
SER_my_4_SR_L3 = SER_my_4_SR_L3/sim;

figure
semilogy(power_range,SER_my_SR_L1,'--o','color',[0.8,0,0],'LineWidth',1.5)
hold on
semilogy(power_range,SER_my_2_SR_L1,'--s','color',[0,0.4,0.8],'LineWidth',1.5)
semilogy(power_range,SER_my_4_SR_L1,'--+','color',[0.15,0.15,0.15],'LineWidth',1.5)
semilogy(power_range,SER_my_SR_L2,'-o','color',[0.8,0,0],'LineWidth',1.5)
semilogy(power_range,SER_my_2_SR_L2,'-s','color',[0,0.4,0.8],'LineWidth',1.5)
semilogy(power_range,SER_my_4_SR_L2,'-+','color',[0.15,0.15,0.15],'LineWidth',1.5)
semilogy(power_range,SER_my_SR_L3,'-.o','color',[0.8,0,0],'LineWidth',1.5)
semilogy(power_range,SER_my_2_SR_L3,'-.s','color',[0,0.4,0.8],'LineWidth',1.5)
semilogy(power_range,SER_my_4_SR_L3,'-.+','color',[0.15,0.15,0.15],'LineWidth',1.5)
hold off
xlabel('Transmit power {\itP} (dBm)');
ylabel('SRs SER');
grid on
axis([22 28 1e-6 0.05])
legend('Continuous, {\it L}=8','2-bit, {\it L}=8','4-bit, {\it L}=8',...
    'Continuous, {\it L}=12','2-bit, {\it L}=12','4-bit, {\it L}=12',...
    'Continuous, {\it L}=16','2-bit, {\it L}=16','4-bit, {\it L}=16')


