% This Matlab script can be used to generate Fig. 4 and Fig. 5 in the paper:
% R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

% cvx_solver mosek

clear;
clc;

load('Channel.mat');

P_dBm_range = (8:1:14);
alpha_range = (2:0.2:3.2);
%%%  average SER
SER_CE_a = zeros(N_sim,length(P_dBm_range)); %%% continuous
SER_CE_q2_a = zeros(N_sim,length(P_dBm_range)); %%% quantize, 2-bit
SER_CE_q3_a = zeros(N_sim,length(P_dBm_range)); %%% quantize, 3-bit
SER_CE_q4_a = zeros(N_sim,length(P_dBm_range)); %%% quantize, 4-bit
SER_CE_q5_a = zeros(N_sim,length(P_dBm_range)); %%% quantize, 5-bit
SER_CE_g2_a = zeros(N_sim,length(P_dBm_range)); %%% heuristic, 2-bit
SER_CE_g3_a = zeros(N_sim,length(P_dBm_range)); %%% heuristic, 3-bit
SER_CE_g4_a = zeros(N_sim,length(P_dBm_range)); %%% heuristic, 4-bit
SER_CE_g5_a = zeros(N_sim,length(P_dBm_range)); %%% heuristic, 5-bit
SER_CE_o1_a = zeros(N_sim,length(P_dBm_range)); %%% BnB, 1bit
SER_CE_o2_a = zeros(N_sim,length(P_dBm_range)); %%% BnB, 2bit
%%% maximize SER
SER_CE_m = zeros(N_sim,length(P_dBm_range));
SER_CE_q2_m = zeros(N_sim,length(P_dBm_range));
SER_CE_q3_m = zeros(N_sim,length(P_dBm_range));
SER_CE_q4_m = zeros(N_sim,length(P_dBm_range));
SER_CE_q5_m = zeros(N_sim,length(P_dBm_range));
SER_CE_g2_m = zeros(N_sim,length(P_dBm_range));
SER_CE_g3_m = zeros(N_sim,length(P_dBm_range));
SER_CE_g4_m = zeros(N_sim,length(P_dBm_range));
SER_CE_g5_m = zeros(N_sim,length(P_dBm_range));
SER_CE_o1_m = zeros(N_sim,length(P_dBm_range));
SER_CE_o2_m = zeros(N_sim,length(P_dBm_range));
%%% average transmit power
power_CE = zeros(N_sim,length(alpha_range));
power_CE_q2 = zeros(N_sim,length(alpha_range));
power_CE_q3 = zeros(N_sim,length(alpha_range));
power_CE_q4 = zeros(N_sim,length(alpha_range));
power_CE_q5 = zeros(N_sim,length(alpha_range));
power_CE_g2 = zeros(N_sim,length(alpha_range));
power_CE_g3 = zeros(N_sim,length(alpha_range));
power_CE_g4 = zeros(N_sim,length(alpha_range));
power_CE_g5 = zeros(N_sim,length(alpha_range));
power_CE_o1 = zeros(N_sim,length(alpha_range));
power_CE_o2 = zeros(N_sim,length(alpha_range));

for sim = 1:N_sim
    sim
    tic
    H = H_tot(:,:,sim);

    Noise = sqrt(sigma2/2)*(randn(K,Ns)+1i*randn(K,Ns));
    S = exp(1i*2*pi*rand(K,Ns));
    [S,index_S] = get_adaptive_modulate(S,4*ones(K,1));

    [X_CE,t] = get_X_CE(H);

    %%%% quantize
    X_CE_q2 = exp( 1i*delta2.*round( angle(X_CE)./delta2 ) );
    t_q2 = min(min((real(H*X_CE_q2.*exp(-1i.*phi_u))-abs(imag(H*X_CE_q2.*exp(-1i.*phi_u))))/sqrt(2)));
    X_CE_q3 = exp( 1i*delta3.*round( angle(X_CE)./delta3 ) );
    t_q3 = min(min((real(H*X_CE_q3.*exp(-1i.*phi_u))-abs(imag(H*X_CE_q3.*exp(-1i.*phi_u))))/sqrt(2)));
    X_CE_q4 = exp( 1i*delta4.*round( angle(X_CE)./delta4 ) );
    t_q4 = min(min((real(H*X_CE_q4.*exp(-1i.*phi_u))-abs(imag(H*X_CE_q4.*exp(-1i.*phi_u))))/sqrt(2)));
    X_CE_q5 = exp( 1i*delta5.*round( angle(X_CE)./delta5 ) );
    t_q5 = min(min((real(H*X_CE_q5.*exp(-1i.*phi_u))-abs(imag(H*X_CE_q5.*exp(-1i.*phi_u))))/sqrt(2)));

    %%%% greedy
    [X_CE_g2,t_g2] = get_X_greedy(H,X_CE,B2);
    [X_CE_g3,t_g3] = get_X_greedy(H,X_CE,B3);
    [X_CE_g4,t_g4] = get_X_greedy(H,X_CE,B4);
    [X_CE_g5,t_g5] = get_X_greedy(H,X_CE,B5);

    %%%% branch-and-bound
    [X_CE_o1,t_o1] = get_X_cvx(H,B1);
    [X_CE_o2,t_o2] = get_X_cvx(H,B2);

    power_CE(sim,:) = 30 - 20*log10(abs(t)) + 20.*log10(alpha_range) + 10*log10(sigma2);
    power_CE_q2(sim,:) = 30 - 20*log10(abs(t_q2)) + 20.*log10(alpha_range) + 10*log10(sigma2);
    power_CE_q3(sim,:) = 30 - 20*log10(abs(t_q3)) + 20.*log10(alpha_range) + 10*log10(sigma2);
    power_CE_q4(sim,:) = 30 - 20*log10(abs(t_q4)) + 20.*log10(alpha_range) + 10*log10(sigma2);
    power_CE_q5(sim,:) = 30 - 20*log10(abs(t_q5)) + 20.*log10(alpha_range) + 10*log10(sigma2);
    power_CE_g2(sim,:) = 30 - 20*log10(abs(t_g2)) + 20.*log10(alpha_range) + 10*log10(sigma2);
    power_CE_g3(sim,:) = 30 - 20*log10(abs(t_g3)) + 20.*log10(alpha_range) + 10*log10(sigma2);
    power_CE_g4(sim,:) = 30 - 20*log10(abs(t_g4)) + 20.*log10(alpha_range) + 10*log10(sigma2);
    power_CE_g5(sim,:) = 30 - 20*log10(abs(t_g5)) + 20.*log10(alpha_range) + 10*log10(sigma2);
    power_CE_o1(sim,:) = 30 - 20*log10(abs(t_o1)) + 20.*log10(alpha_range) + 10*log10(sigma2);
    power_CE_o2(sim,:) = 30 - 20*log10(abs(t_o2)) + 20.*log10(alpha_range) + 10*log10(sigma2);

    for P_index = 1:length(P_dBm_range)
        P = P_dBm_range(P_index);

        R = sqrt(10^(0.1*P-3))*H*X_CE(:,index_S) + Noise;
        R_q2 = sqrt(10^(0.1*P-3))*H*X_CE_q2(:,index_S) + Noise;
        R_q3 = sqrt(10^(0.1*P-3))*H*X_CE_q3(:,index_S) + Noise;
        R_q4 = sqrt(10^(0.1*P-3))*H*X_CE_q4(:,index_S) + Noise;
        R_q5 = sqrt(10^(0.1*P-3))*H*X_CE_q5(:,index_S) + Noise;
        R_g2 = sqrt(10^(0.1*P-3))*H*X_CE_g2(:,index_S) + Noise;
        R_g3 = sqrt(10^(0.1*P-3))*H*X_CE_g3(:,index_S) + Noise;
        R_g4 = sqrt(10^(0.1*P-3))*H*X_CE_g4(:,index_S) + Noise;
        R_g5 = sqrt(10^(0.1*P-3))*H*X_CE_g5(:,index_S) + Noise;
        R_o1 = sqrt(10^(0.1*P-3))*H*X_CE_o1(:,index_S) + Noise;
        R_o2 = sqrt(10^(0.1*P-3))*H*X_CE_o2(:,index_S) + Noise;

        temp1 = get_SER(R,S,4.*ones(K,1));
        SER_CE_a(sim,P_index) = sum(temp1)/K;
        SER_CE_m(sim,P_index) = max(temp1);
        temp3 = get_SER(R_q2,S,4.*ones(K,1));
        SER_CE_q2_a(sim,P_index) = sum(temp3)/K;
        SER_CE_q2_m(sim,P_index) = max(temp3);
        temp4 = get_SER(R_q3,S,4.*ones(K,1));
        SER_CE_q3_a(sim,P_index) = sum(temp4)/K;
        SER_CE_q3_m(sim,P_index) = max(temp4);
        temp5 = get_SER(R_q4,S,4.*ones(K,1));
        SER_CE_q4_a(sim,P_index) = sum(temp5)/K;
        SER_CE_q4_m(sim,P_index) = max(temp5);
        temp6 = get_SER(R_q5,S,4.*ones(K,1));
        SER_CE_q5_a(sim,P_index) = sum(temp6)/K;
        SER_CE_q5_m(sim,P_index) = max(temp6);
        temp8 = get_SER(R_g2,S,4.*ones(K,1));
        SER_CE_g2_a(sim,P_index) = sum(temp8)/K;
        SER_CE_g2_m(sim,P_index) = max(temp8);
        temp9 = get_SER(R_g3,S,4.*ones(K,1));
        SER_CE_g3_a(sim,P_index) = sum(temp9)/K;
        SER_CE_g3_m(sim,P_index) = max(temp9);
        temp10 = get_SER(R_g4,S,4.*ones(K,1));
        SER_CE_g4_a(sim,P_index) = sum(temp10)/K;
        SER_CE_g4_m(sim,P_index) = max(temp10);
        temp11 = get_SER(R_g5,S,4.*ones(K,1));
        SER_CE_g5_a(sim,P_index) = sum(temp11)/K;
        SER_CE_g5_m(sim,P_index) = max(temp11);
        temp12 = get_SER(R_o1,S,4.*ones(K,1));
        SER_CE_o1_a(sim,P_index) = sum(temp12)/K;
        SER_CE_o1_m(sim,P_index) = max(temp12);
        temp13 = get_SER(R_o2,S,4.*ones(K,1));
        SER_CE_o2_a(sim,P_index) = sum(temp13)/K;
        SER_CE_o2_m(sim,P_index) = max(temp13);
    end
    toc
end

SER_CE_a = SER_CE_a/sim;
SER_CE_q2_a = SER_CE_q2_a/sim;
SER_CE_q3_a = SER_CE_q3_a/sim;
SER_CE_q4_a = SER_CE_q4_a/sim;
SER_CE_q5_a = SER_CE_q5_a/sim;
SER_CE_g2_a = SER_CE_g2_a/sim;
SER_CE_g3_a = SER_CE_g3_a/sim;
SER_CE_g4_a = SER_CE_g4_a/sim;
SER_CE_g5_a = SER_CE_g5_a/sim;
SER_CE_o1_a = SER_CE_o1_a/sim;
SER_CE_o2_a = SER_CE_o2_a/sim;
SER_CE_m = SER_CE_m/sim;
SER_CE_q2_m = SER_CE_q2_m/sim;
SER_CE_q3_m = SER_CE_q3_m/sim;
SER_CE_q4_m = SER_CE_q4_m/sim;
SER_CE_q5_m = SER_CE_q5_m/sim;
SER_CE_g2_m = SER_CE_g2_m/sim;
SER_CE_g3_m = SER_CE_g3_m/sim;
SER_CE_g4_m = SER_CE_g4_m/sim;
SER_CE_g5_m = SER_CE_g5_m/sim;
SER_CE_o1_m = SER_CE_o1_m/sim;
SER_CE_o2_m = SER_CE_o2_m/sim;
power_CE = power_CE/sim;
power_CE_q2 = power_CE_q2/sim;
power_CE_q3 = power_CE_q3/sim;
power_CE_q4 = power_CE_q4/sim;
power_CE_q5 = power_CE_q5/sim;
power_CE_g2 = power_CE_g2/sim;
power_CE_g3 = power_CE_g3/sim;
power_CE_g4 = power_CE_g4/sim;
power_CE_g5 = power_CE_g5/sim;
power_CE_o1 = power_CE_o1/sim;
power_CE_o2 = power_CE_o2/sim;

figure
plot(alpha_range,power_CE,'-o','color',[0.8,0,0],'LineWidth',1.5)
hold on
plot(alpha_range,power_CE_q2,'-s','color',[0,0.4,0.8],'LineWidth',1.5)
plot(alpha_range,power_CE_q3,'-d','color',[0,0.7,0.7],'LineWidth',1.5)
plot(alpha_range,power_CE_q4,'-+','color',[0.15,0.15,0.15],'LineWidth',1.5)
plot(alpha_range,power_CE_q5,'-x','color',[0,0.5,0],'LineWidth',1.5)
plot(alpha_range,power_CE_g2,'--s','color',[0,0.4,0.8],'LineWidth',1.5)
plot(alpha_range,power_CE_g3,'--d','color',[0,0.7,0.7],'LineWidth',1.5)
plot(alpha_range,power_CE_g4,'--+','color',[0.15,0.15,0.15],'LineWidth',1.5)
plot(alpha_range,power_CE_g5,'--x','color',[0,0.5,0],'LineWidth',1.5)
plot(alpha_range,power_CE_o1,':','color',[0.7,0.2,0.6],'LineWidth',1.5)
plot(alpha_range,power_CE_o2,':','color',[0,0.4,0.8],'LineWidth',1.5)
hold off
xlabel('QoS requirement \alpha');
ylabel('Average transmit power (dBm)');
grid on
set(gca,'XTicklabel',{'2\sigma','2.2\sigma','2.4\sigma','2.6\sigma','2.8\sigma','3\sigma','3.2\sigma',})
legend('Continuous','Quantize, 2-bit','Quantize, 3-bit','Quantize, 4-bit',...
    'Quantize, 5-bit','Heuristic, 2-bit','Heuristic, 3-bit','Heuristic, 4-bit',...
    'Heuristic, 5-bit','B & B, 1-bit','B & B, 2-bit');%,'o2');
axis([2 3.2 9.5 20])

figure
semilogy(P_dBm_range,SER_CE_a,'-o','color',[0.8,0,0],'LineWidth',1.5)
hold on
semilogy(P_dBm_range,SER_CE_o2_a,':','color',[0,0.4,0.8],'LineWidth',1.5)
semilogy(P_dBm_range,SER_CE_g3_a,'--d','color',[0,0.7,0.7],'LineWidth',1.5)
semilogy(P_dBm_range,SER_CE_g4_a,'--+','color',[0.15,0.15,0.15],'LineWidth',1.5)
semilogy(P_dBm_range,SER_CE_q5_a,'-x','color',[0,0.5,0],'LineWidth',1.5)
hold off
xlabel('Average transmit power {\it P} (dBm)');
ylabel('Average SER');
grid on
axis([8 14 1e-6 0.03])
legend('Continuous','B & B, 2-bit','Heuristic, 3-bit','Heuristic, 4-bit',...
    'Quantize, 5-bit');

figure
semilogy(P_dBm_range,SER_CE_m,'-o','color',[0.8,0,0],'LineWidth',1.5)
hold on
semilogy(P_dBm_range,SER_CE_o2_m,':','color',[0,0.4,0.8],'LineWidth',1.5)
semilogy(P_dBm_range,SER_CE_g3_m,'--d','color',[0,0.7,0.7],'LineWidth',1.5)
semilogy(P_dBm_range,SER_CE_g4_m,'--+','color',[0.15,0.15,0.15],'LineWidth',1.5)
semilogy(P_dBm_range,SER_CE_q5_m,'-x','color',[0,0.5,0],'LineWidth',1.5)
hold off
xlabel('Average transmit power {\it P} (dBm)');
ylabel('Maximum SER');
grid on
axis([8 14 1e-6 0.03])
legend('Continuous','B & B, 2-bit','Heuristic, 3-bit','Heuristic, 4-bit',...
    'Quantize, 5-bit');
