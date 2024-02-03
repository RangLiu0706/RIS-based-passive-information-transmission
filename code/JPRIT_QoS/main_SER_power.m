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
power_range = (20:1:26);

SER_my_PR = zeros(1,length(power_range));
SER_my_2_PR = zeros(1,length(power_range));
SER_my_3_PR = zeros(1,length(power_range));
SER_my_4_PR = zeros(1,length(power_range));
SER_my_5_PR = zeros(1,length(power_range));
SER_wo_SR = zeros(1,length(power_range));
SER_BS_PR = zeros(1,length(power_range));

B2 = 2;
B3 = 3;
B4 = 4;
B5 = 5;
d_ar = 10;  %% BS-IRS distance
d_rs = 20; %%%%
d_ru = 100;

belta1 = sqrt(10^(0.3)/(1+10^(0.3)));
belta2 = sqrt(1/(1+10^(0.3)));

H_au = zeros(K+1,M);
H_ru = zeros(K+1,N);

Nmax = 9;
res_th = 1e-3;
Ns = 100000;
L = 10;
R = zeros(K+1,Ns);
R2 = zeros(K+1,Ns);
R3 = zeros(K+1,Ns);
R4 = zeros(K+1,Ns);
R5 = zeros(K+1,Ns);
for sim = 1:N_sim
    sim
    tic
    d_as = sqrt(d_ar^2+d_rs^2)+(d_ar+d_rs-sqrt(d_ar^2+d_rs^2))*rand;
    d_au = sqrt(d_ar^2+d_ru^2)+(d_ar+d_ru-sqrt(d_ar^2+d_ru^2)).*rand(1,K);

    H_ar = sqrt(10^(-3)*d_ar^(-2))*(belta1*channel_ar(M,N)+belta2*(randn(N,M)+1i*randn(N,M))/sqrt(2));
    H_ru(1:K,:) = sqrt(10^(-3)*d_ru^(-3))*(randn(K,N)+1i*randn(K,N))/sqrt(2);
    H_ru(K+1,:) = sqrt(10^(-3)*d_rs^(-3))*(randn(1,N)+1i*randn(1,N))/sqrt(2);
    for i = 1:1:K
        H_au(i,:) = sqrt(10^(-3)*d_au(i)^(-3))*(randn(1,M)+1i*randn(1,M))/sqrt(2);
    end
    H_au(K+1,:) = sqrt(10^(-3)*d_as^(-3))*(randn(1,M)+1i*randn(1,M))/sqrt(2);

    Noise = sqrt(0.5*sigma2)*(randn(K+1,Ns)+1i*randn(K+1,Ns));
    S_PR = exp(1i*2*pi*rand(K,Ns));
    S_SR = exp(1i*2*pi*rand(1,Ns/L));
    temp = reshape(repmat(S_SR.',1,L).',Ns,1);
    [S_BS,index_BS] = get_adaptive_modulate([S_PR;temp.'],[omega 2]);
    [S_PR,index_PR] = get_adaptive_modulate(S_PR,omega);  %% modulate PRs' information
    [S_SR,index_SR] = get_adaptive_modulate(S_SR,2);  %% modulate SR's information
    index_theta = reshape(repmat(index_SR,L,1),Ns,1)';

    for power_index = 1:length(power_range)
        power = power_range(power_index)

        [X_my,Theta_my,t] = get_X_theta(H_au,H_ar,H_ru,dist,power,Nmax,res_th);
        [X_my2,Theta_my2,t2] = get_X_theta_b(H_au,H_ar,H_ru,dist,power,B2,Nmax,res_th);
        [X_my3,Theta_my3,t3] = get_X_theta_b(H_au,H_ar,H_ru,dist,power,B3,Nmax,res_th);
        [X_my4,Theta_my4,t4] = get_X_theta_b(H_au,H_ar,H_ru,dist,power,B4,Nmax,res_th);
        [X_my5,Theta_my5,t5] = get_X_theta_b(H_au,H_ar,H_ru,dist,power,B5,Nmax,res_th);
        [X6,theta6,t6] = get_X_theta_wo_SR(H_au(1:K,:),H_ar,H_ru(1:K,:),dist(1:K),omega,power,Nmax,res_th);
        [X7,theta7,t7] = get_X_theta_wo_SR(H_au,H_ar,H_ru,dist,[omega 2],power,Nmax,res_th);

        for i = 1:1:Ns
            R(:,i) = (H_au + H_ru*diag(Theta_my(:,index_theta(i)))*H_ar)*X_my(:,index_PR(i))+ Noise(:,i);
            R2(:,i) = (H_au + H_ru*diag(Theta_my2(:,index_theta(i)))*H_ar)*X_my2(:,index_PR(i))+ Noise(:,i);
            R3(:,i) = (H_au + H_ru*diag(Theta_my3(:,index_theta(i)))*H_ar)*X_my3(:,index_PR(i))+ Noise(:,i);
            R4(:,i) = (H_au + H_ru*diag(Theta_my4(:,index_theta(i)))*H_ar)*X_my4(:,index_PR(i))+ Noise(:,i);
            R5(:,i) = (H_au + H_ru*diag(Theta_my5(:,index_theta(i)))*H_ar)*X_my5(:,index_PR(i))+ Noise(:,i);
        end
        SER_my_PR(power_index) = SER_my_PR(power_index) + sum(get_SER(R(1:K,:),S_PR,omega))/K;
        SER_my_2_PR(power_index) = SER_my_2_PR(power_index) + sum(get_SER(R2(1:K,:),S_PR,omega))/K;
        SER_my_3_PR(power_index) = SER_my_3_PR(power_index) + sum(get_SER(R3(1:K,:),S_PR,omega))/K;
        SER_my_4_PR(power_index) = SER_my_4_PR(power_index) + sum(get_SER(R4(1:K,:),S_PR,omega))/K;
        SER_my_5_PR(power_index) = SER_my_5_PR(power_index) + sum(get_SER(R5(1:K,:),S_PR,omega))/K;

        R6 = ( H_au(1:K,:) + H_ru(1:K,:)*diag(theta6)*H_ar )*X6(:,index_PR) + Noise(1:K,:);
        R7 = ( H_au + H_ru*diag(theta7)*H_ar )*X7(:,index_BS) + Noise;
        SER_wo_SR(power_index) = SER_wo_SR(power_index) + sum(get_SER(R6,S_PR,omega))/K;
        SER_BS_PR(power_index) = SER_BS_PR(power_index) + sum(get_SER(R7(1:K,:),S_BS(1:K,:),omega))/K;

    end
    toc
end

SER_my_PR = SER_my_PR/sim;
SER_my_2_PR = SER_my_2_PR/sim;
SER_my_3_PR = SER_my_3_PR/sim;
SER_my_4_PR = SER_my_4_PR/sim;
SER_my_5_PR = SER_my_5_PR/sim;
SER_BS_PR = SER_BS_PR/N_sim;
SER_wo_SR = SER_wo_SR/sim;

figure
semilogy(power_range,SER_my_PR,'-o','color',[0.85,0.1,0.1],'LineWidth',1.5)
hold on
semilogy(power_range,SER_my_2_PR,'-d','color',[0.1,0.1,0.1],'LineWidth',1.5)
semilogy(power_range,SER_my_3_PR,'-s','color',[0.54,0.3,0.35],'LineWidth',1.5)
semilogy(power_range,SER_my_4_PR,'-^','color',[0,0.5,0],'LineWidth',1.5)
semilogy(power_range,SER_my_5_PR,'-^','color',[0,0.5,0],'LineWidth',1.5)
semilogy(power_range,SER_wo_SR,'-+','color',[0.8,0.5,0],'LineWidth',1.5,'markersize',7)
semilogy(power_range,SER_BS_PR,'-x','color',[0.1,0.5,0.6],'LineWidth',1.5,'markersize',7)
hold off
xlabel('Transmit power {\itP} (dBm)');
ylabel('PRs average SER');
grid on
% axis([15 20 3.0e-5 0.1])
legend('Continuous','2-bit','3-bit','4-bit','5-bit','w/o SR','BS, w/ SR');

