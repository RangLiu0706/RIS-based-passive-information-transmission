% Solve for the problem (45).
% This is used in the paper: R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

function [X,Theta,tt] = get_X_theta(H_au,H_ar,H_ru,dist,power,Nmax,res_th)

global sigma2
[K,M] = size(H_au);
K = K - 1;
[~,N] = size(H_ru);
a = zeros(2*(K+1)*4^K,1);
B = zeros(2*(K+1)*4^K,2*N);
phi_u = zeros(K,4^K);
t1 = [1;0];
t2 = [0;1];

for k = 0:1:4^K-1
    s = dec2bin(k,2*K);
    for i = 1:1:K
        temp = 2*(s(2*i-1)-48)+s(2*i)-48;
        phi_u(i,k+1) = pi/4 + temp*pi/2;
    end
end

%%% random initialize
temp = exp(1i*2*pi*rand(N,1));
theta = [temp;conj(temp)];
Theta = [theta(1:N) theta(N+1:end)];
H1 = (H_au + H_ru*diag(Theta(:,1))*H_ar)./(repmat(dist.'*sqrt(sigma2),1,M));
H2 = (H_au + H_ru*diag(Theta(:,2))*H_ar)./(repmat(dist.'*sqrt(sigma2),1,M));

cvx_begin quiet
variable X(M,4^K) complex
variable t
maximize t
subject to
r1 = H1(1:K,:)*X.*exp(-1i*phi_u);
r2 = H2(1:K,:)*X.*exp(-1i*phi_u);
t <= real(r1)*sin(pi/4) - abs(imag(r1))*cos(pi/4);
t <= real(r2)*sin(pi/4) - abs(imag(r2))*cos(pi/4);
t <= imag(H1(K+1,:)*X);
t <= -imag(H2(K+1,:)*X);
square_pos(norm(X,'fro')) <= 4^K*10^(0.1*power-3);
cvx_end

res = 1;
iter = 1;
tt(1,iter) = t;
while iter <= Nmax && res >= res_th
    t_re = t;

    for k = 1:1:K
        for j = 1:1:4^K
            a(2*(k-1)*4^K+2*j-1,1) = H_au(k,:)*X(:,j)*exp(-1i*phi_u(k,j))./dist(k);
            a(2*(k-1)*4^K+2*j,1) = a(2*(k-1)*4^K+2*j-1,1);
            B(2*(k-1)*4^K+2*j-1,:) = kron( t1',H_ru(k,:)*diag(H_ar*X(:,j))*exp(-1i*phi_u(k,j)) )./dist(k);
            B(2*(k-1)*4^K+2*j,:) = kron(t2',H_ru(k,:)*diag(H_ar*X(:,j))*exp(-1i*phi_u(k,j)))./dist(k);
        end
    end
    for j = 1:1:4^K
        a(2*K*4^K+2*j-1,1) = H_au(K+1,:)*X(:,j)*exp(-1i*pi/2)./dist(end);
        a(2*K*4^K+2*j,1) = H_au(K+1,:)*X(:,j)*exp(-1i*3*pi/2)./dist(end);
        B(2*K*4^K+2*j-1,:) = kron(t1',H_ru(K+1,:)*diag(H_ar*X(:,j))*exp(-1i*pi/2))./dist(end);
        B(2*K*4^K+2*j,:) = kron(t2',H_ru(K+1,:)*diag(H_ar*X(:,j))*exp(-1i*3*pi/2))./dist(end);
    end

    [theta] = get_theta(a.*1e5,B.*1e5,K,theta);
    Theta = [theta(1:N) theta(N+1:2*N)];

    H1 = (H_au + H_ru*diag(Theta(:,1))*H_ar)./(repmat(dist.'*sqrt(sigma2),1,M));
    H2 = (H_au + H_ru*diag(Theta(:,2))*H_ar)./(repmat(dist.'*sqrt(sigma2),1,M));

    cvx_begin quiet
    variable X(M,4^K) complex
    variable t
    maximize t
    subject to
    r1 = H1(1:K,:)*X.*exp(-1i*phi_u);
    r2 = H2(1:K,:)*X.*exp(-1i*phi_u);
    t <= real(r1)*sin(pi/4) - abs(imag(r1))*cos(pi/4);
    t <= real(r2)*sin(pi/4) - abs(imag(r2))*cos(pi/4);
    t <= imag(H1(K+1,:)*X);
    t <= -imag(H2(K+1,:)*X);
    square_pos(norm(X,'fro')) <= 4^K*10^(0.1*power-3);
    cvx_end

    res = abs(1-t/t_re);
    iter = iter + 1;
    tt(1,iter) = t;
end

end