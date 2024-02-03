% The proposed heuristic algorithm to obtain low-resolution results.
% This is used in the paper: R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02
% Inputs: H: the channel; X_CE: the continuous results; B: the resolution;
% Outputs: X_g: RIS reflecting coefficients;
%          t_g: the objective value

function [X_g,t_g] = get_X_greedy(H,X_CE,B)
[K,N] = size(H);
delta = 2*pi/(2^B);
X_g = X_CE;
phi_u = zeros(K,4^K);
A = zeros(K,N);
x_set = exp(1i*delta*(1:1:2^B));

Nmax = 100;
res_th = 1e-3;

for m = 0:1:4^K-1
    s = dec2bin(m,2*K);
    for i = 1:1:K
        temp = 2*(s(2*i-1)-48)+s(2*i)-48;
        phi_u(i,m+1) = pi/4 + temp*pi/2;
        A(i,:) = H(i,:).*exp(-1i*phi_u(i,m+1));
    end
    
    xm = X_g(:,m+1);
    iter = 1;
    res = 1;
    obj = min(real(A*xm)-abs(imag(A*xm)));
    while iter <= Nmax && res >= res_th
        objpre = obj;
        for i = 1:1:N
            xm(i) = 0;
            temp = repmat(A*xm,1,2^B) + repmat(A(:,i),1,2^B).*repmat(x_set,K,1);
            dec = min(real(temp)-abs(imag(temp)));
            [~,ind] = find(dec == max(dec),1);
            xm(i) = x_set(ind);
        end
        obj = min(real(A*xm)-abs(imag(A*xm)));
        res = abs(1-objpre/obj);
        iter = iter + 1;
    end
    X_g(:,m+1) = xm;    
end
t_g = min(min((real(H*X_g.*exp(-1i.*phi_u))-abs(imag(H*X_g.*exp(-1i.*phi_u))))/sqrt(2)));



