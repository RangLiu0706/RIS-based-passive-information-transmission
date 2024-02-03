% This Matlab script can be used to generate Fig. 9 in the paper:
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

dist_SR = 0.5;  %% QoS requirement of secondary receiver
dist_PR = 2.5;  %% QoS of PRs
dist = [dist_PR.*ones(1,K) dist_SR];
omega = 4.*ones(1,K);
power = 25;  %% dBm

N_sim = 1000;  %% number of simulations
iter_range = (1:1:20);

t_my = zeros(1,length(iter_range));
t_my_1 = zeros(1,length(iter_range));
t_my_2 = zeros(1,length(iter_range));
t_my_3 = zeros(1,length(iter_range));
t_my_4 = zeros(1,length(iter_range));
t_my_5 = zeros(1,length(iter_range));

B1 = 1;
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

    [X_my,Theta_my,t] = get_X_theta(H_au,H_ar,H_ru,dist,power,Nmax,res_th);
    [X_my1,Theta_my1,t1] = get_X_theta_b(H_au,H_ar,H_ru,dist,power,B1,Nmax,res_th);
    [X_my2,Theta_my2,t2] = get_X_theta_b(H_au,H_ar,H_ru,dist,power,B2,Nmax,res_th);
    [X_my3,Theta_my3,t3] = get_X_theta_b(H_au,H_ar,H_ru,dist,power,B3,Nmax,res_th);
    [X_my4,Theta_my4,t4] = get_X_theta_b(H_au,H_ar,H_ru,dist,power,B4,Nmax,res_th);
    [X_my5,Theta_my5,t5] = get_X_theta_b(H_au,H_ar,H_ru,dist,power,B5,Nmax,res_th);

    t_my = t_my + t;
    t_my_1 = t_my_1 + t1;
    t_my_2 = t_my_2 + t2;
    t_my_3 = t_my_3 + t3;
    t_my_4 = t_my_4 + t4;
    t_my_5 = t_my_5 + t5;
    toc
end

t_my = t_my/sim;
t_my_1 = t_my_1/sim;
t_my_2 = t_my_2/sim;
t_my_3 = t_my_3/sim;
t_my_4 = t_my_4/sim;
t_my_5 = t_my_5/sim;

figure
plot(iter_range,t_my,'-o','color',[0.5,0,0],'LineWidth',1.5)
hold on
plot(iter_range,t_my_1,'->','color',[0.5,0.5,0],'LineWidth',1.5)
plot(iter_range,t_my_2,'-d','color',[0,0.5,0],'LineWidth',1.5)
plot(iter_range,t_my_3,'-^','color',[0,0,0.5],'LineWidth',1.5)
plot(iter_range,t_my_4,'-s','color',[0,0.5,0.5],'LineWidth',1.5)
plot(iter_range,t_my_5,'-+','color',[0.1,0.1,0.1],'LineWidth',1.5)
hold off
xlabel('Number of iterations');
ylabel('Average transmit power (dBm)');
grid on
legend('Proposed, continuous','Proposed, 1-bit','Proposed, 2-bit',...
    'Proposed, 3-bit','Proposed, 4-bit','Proposed, 5-bit')
