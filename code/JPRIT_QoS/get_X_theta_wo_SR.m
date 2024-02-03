% Comparison scheme.
% This is used in the paper: R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

function [X,theta,tt] = get_X_theta_wo_SR(H_au,H_ar,H_ru,dist,omega,power,Nmax,res_th)
global sigma2;
[K,M] = size(H_au);
[~,N] = size(H_ru);
nt = prod(omega);
a = zeros(K*nt,1);
B = zeros(K*nt,N);
phi_u = zeros(K,nt);
phi = pi./omega';
Phi = repmat(phi,1,nt);
Dist = sqrt(sigma2).*repmat(dist',1,nt);
if omega(end) == 2
    ws = 1;
else
    ws = 0;
end
for m = 0:1:nt-1
    s = dec2bin(m,log2(nt));
    indicator = 1;
    for i = 1:1:K
        temp = 0;
        for j = 1:1:log2(omega(i))
            temp = temp + 2^(log2(omega(i))-j)*(s(indicator)-48);
            indicator = indicator + 1;
        end
        phi_u(i,m+1) = pi/omega(i)+temp*2*pi/(omega(i));
    end
end

theta = exp(1i*2*pi*rand(N,1));
Theta = diag(theta);
H = H_au+H_ru*Theta*H_ar;
cvx_begin quiet
variable X(M,nt) complex
variable t
maximize t
subject to
r = H*X.*exp(-1i*phi_u)./Dist;
t <= real(r).*sin(Phi) - abs(imag(r)).*cos(Phi);
square_pos(norm(X,'fro')) <= nt*10^(0.1*power-3);
cvx_end

res = 1;
iter = 1;
tt(1,1) = t;
while iter <= Nmax && res >= res_th
    t_re = t;

    for k = 1:1:K
        for j = 1:1:nt
            a(nt*(k-1)+j) = H_au(k,:)*X(:,j)*exp(-1i*phi_u(k,j))./dist(k);
            B(nt*(k-1)+j,:) = H_ru(k,:)*diag(H_ar*X(:,j)*exp(-1i*phi_u(k,j)))./dist(k);
        end
    end

    theta = get_theta_ws(a.*1e5,B.*1e5,K-1,theta,ws);
    Theta = diag(theta);

    H = H_au+H_ru*Theta*H_ar;
    cvx_begin quiet
    variable X(M,nt) complex
    variable t
    maximize t
    subject to
    r = H*X.*exp(-1i*phi_u)./Dist;
    t <= real(r).*sin(Phi) - abs(imag(r)).*cos(Phi);
    square_pos(norm(X,'fro')) <= nt*10^(0.1*power-3);
    cvx_end

    res = abs(1-t/t_re);
    iter = iter + 1;
    tt(1,iter) = t;

end
end