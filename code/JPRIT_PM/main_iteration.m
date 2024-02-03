% This Matlab script can be used to generate Fig. 6 in the paper:
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

dist_SR = 0.5;
dist = 2.5.*ones(1,K);
omega = 4.*ones(1,K);

N_sim = 1000;
iter_range = (1:1:20);

power_my = zeros(1,length(iter_range));
power_my_1 = zeros(1,length(iter_range));
power_my_2 = zeros(1,length(iter_range));
power_my_3 = zeros(1,length(iter_range));
power_my_4 = zeros(1,length(iter_range));
power_my_5 = zeros(1,length(iter_range));

B1 = 1;
B2 = 2;
B3 = 3;
B4 = 4;
B5 = 5;

d_ar = 10;
d_rs = 20;
d_ru = 100;
belta1 = sqrt(10^(0.3)/(1+10^(0.3)));
belta2 = sqrt(1/(1+10^(0.3)));
H_au = zeros(K+1,M);
H_ru = zeros(K+1,N);

Nmax = 19;
res_th = -1;
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

    [X,Theta,p] = get_X_theta(H_au,H_ar,H_ru,dist,dist_SR,Nmax,res_th);
    [X_my1,Theta_my1,p1] = get_X_theta_b(H_au,H_ar,H_ru,dist,dist_SR,Nmax,res_th,B1);
    [X_my2,Theta_my2,p2] = get_X_theta_b(H_au,H_ar,H_ru,dist,dist_SR,Nmax,res_th,B2);
    [X_my3,Theta_my3,p3] = get_X_theta_b(H_au,H_ar,H_ru,dist,dist_SR,Nmax,res_th,B3);
    [X_my4,Theta_my4,p4] = get_X_theta_b(H_au,H_ar,H_ru,dist,dist_SR,Nmax,res_th,B4);
    [X_my5,Theta_my5,p5] = get_X_theta_b(H_au,H_ar,H_ru,dist,dist_SR,Nmax,res_th,B5);

    power_my = power_my + 10.*log10(1000.*p./4^K);
    power_my_1 = power_my_1 + 10.*log10(1000.*p1./4^K);
    power_my_2 = power_my_2 + 10.*log10(1000.*p2./4^K);
    power_my_3 = power_my_3 + 10.*log10(1000.*p3./4^K);
    power_my_4 = power_my_4 + 10.*log10(1000.*p4./4^K);
    power_my_5 = power_my_5 + 10.*log10(1000.*p5./4^K);

    toc
end

power_my = power_my/sim;
power_my_1 = power_my_1/sim;
power_my_2 = power_my_2/sim;
power_my_3 = power_my_3/sim;
power_my_4 = power_my_4/sim;
power_my_5 = power_my_5/sim;

figure
plot(iter_range,power_my,'-o','color',[0.5,0,0],'LineWidth',1.5)
hold on
plot(iter_range,power_my_1,'->','color',[0.5,0.5,0],'LineWidth',1.5)
plot(iter_range,power_my_2,'-d','color',[0,0.5,0],'LineWidth',1.5)
plot(iter_range,power_my_3,'-^','color',[0,0,0.5],'LineWidth',1.5)
plot(iter_range,power_my_4,'-s','color',[0,0.5,0.5],'LineWidth',1.5)
plot(iter_range,power_my_5,'-+','color',[0.1,0.1,0.1],'LineWidth',1.5)
hold off
xlabel('Number of iterations');
ylabel('Average transmit power (dBm)');
grid on
legend('Continuous','1-bit','2-bit','3-bit','4-bit','5-bit')








