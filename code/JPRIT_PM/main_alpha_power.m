% This Matlab script can be used to generate Fig. 7 in the paper:
% R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02
clear;
clc;

M = 6;
N = 100;
K = 3;
phi = pi/4;
global sigma2
sigma2 = 10^-11;

N_sim = 1000;
dist_PR_range = (2:0.2:3.2);

power_my_2 = zeros(1,length(dist_PR_range));
power_my_3 = zeros(1,length(dist_PR_range));
power_my_4 = zeros(1,length(dist_PR_range));
power_my_5 = zeros(1,length(dist_PR_range));

power_my = zeros(1,length(dist_PR_range));
power_wo_SR = zeros(1,length(dist_PR_range));
power_BS_SR = zeros(1,length(dist_PR_range));

B2 = 2;
B3 = 3;
B4 = 4;
B5 = 5;

dist_SR = 0.5;
omega = 4.*ones(1,K);

d_ar = 10;
d_rs = 20;
d_ru = 100;
belta1 = sqrt(10^(0.3)/(1+10^(0.3)));
belta2 = sqrt(1/(1+10^(0.3)));

H_au = zeros(K+1,M);
H_ru = zeros(K+1,N);

Nmax = 10;
res_th = 1e-3;
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


    for dist_index = 1:length(dist_PR_range)
        dist = dist_PR_range(dist_index).*ones(1,K)

        [X,Theta,p] = get_X_theta(H_au,H_ar,H_ru,dist,dist_SR,Nmax,res_th);

        [X_my2,Theta_my2,p2] = get_X_theta_b(H_au,H_ar,H_ru,dist,dist_SR,Nmax,res_th,B2);
        [X_my3,Theta_my3,p3] = get_X_theta_b(H_au,H_ar,H_ru,dist,dist_SR,Nmax,res_th,B3);
        [X_my4,Theta_my4,p4] = get_X_theta_b(H_au,H_ar,H_ru,dist,dist_SR,Nmax,res_th,B4);
        [X_my5,Theta_my5,p5] = get_X_theta_b(H_au,H_ar,H_ru,dist,dist_SR,Nmax,res_th,B5);
        [X6,Theta6,p6] = get_X_theta_ws_sum(H_au(1:K,:),H_ar,H_ru(1:K,:),dist,omega,Nmax,res_th);
        [X7,Theta7,p7] = get_X_theta_ws_sum(H_au,H_ar,H_ru,[dist dist_SR],[omega 2],Nmax,res_th);


        power_my(dist_index) = power_my(dist_index) + 10*log10(1000*min(p)/4^K);
        power_my_2(dist_index) = power_my_2(dist_index) + 10*log10(1000*min(p2)/4^K);
        power_my_3(dist_index) = power_my_3(dist_index) + 10*log10(1000*min(p3)/4^K);
        power_my_4(dist_index) = power_my_4(dist_index) + 10*log10(1000*min(p4)/4^K);
        power_my_5(dist_index) = power_my_5(dist_index) + 10*log10(1000*min(p5)/4^K);
        power_wo_SR(dist_index) = power_wo_SR(dist_index) + 10*log10(1000*min(p6)/4^K);
        power_BS_SR(dist_index) = power_BS_SR(dist_index) + 10*log10(1000*min(p7)/4^(K+1));
    end
    toc
end

power_my = power_my/sim;
power_my_2 = power_my_2/sim;
power_my_3 = power_my_3/sim;
power_my_4 = power_my_4/sim;
power_my_5 = power_my_5/sim;
power_wo_SR = power_wo_SR/sim;
power_BS_SR = power_BS_SR/sim;

figure
plot(dist_PR_range,power_my,'-o','color',[0.8,0,0],'LineWidth',1.5)
hold on
plot(dist_PR_range,power_my_2,'-s','color',[0,0.4,0.8],'LineWidth',1.5)
plot(dist_PR_range,power_my_3,'-d','color',[0,0.7,0.7],'LineWidth',1.5)
plot(dist_PR_range,power_my_4,'-+','color',[0.15,0.15,0.15],'LineWidth',1.5)
plot(dist_PR_range,power_my_5,'-x','color',[0,0.5,0],'LineWidth',1.5)
plot(dist_PR_range,power_wo_SR,'-^','color',[1,0.3,0.6],'LineWidth',1.5)
plot(dist_PR_range,power_BS_SR-10*log10(0.5),'-*','color',[0.5,0.5,0],'LineWidth',1.5)
hold off
xlabel('PIRs QoS requirement \alpha');
ylabel('Average transmit power (dBm)');
grid on
axis([2 3.2 18.5 25.5])
set(gca,'XTicklabel',{'2\sigma','2.2\sigma','2.4\sigma','2.6\sigma','2.8\sigma','3\sigma','3.2\sigma',})
legend('Continuous','2-bit','3-bit','4-bit', '5-bit','w/o SIR', 'BS, w/ SIR')


