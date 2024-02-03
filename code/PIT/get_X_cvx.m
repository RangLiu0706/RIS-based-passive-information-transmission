% The proposed branch-and-bound algorithm to obtain low-resolution results.
% This is used in the paper: R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02
% Inputs: H: the channel; B: the resolution;
% Outputs: X_o: RIS reflecting coefficients;
%          t_o: the objective value
function [X_o,t_o] = get_X_cvx(H,B)
[K,N] = size(H);
H = H.*10^5.5;
phi_u = zeros(K,4^K);
X_o = zeros(N,4^K);
delta = 2*pi/(2^B);
ind = (1:1:2^B);
v = exp(1i*(delta.*ind'));

for m = 0:1:4^K-1
    
    s = dec2bin(m,2*K);
    for i = 1:1:K
        temp = 2*(s(2*i-1)-48)+s(2*i)-48;
        phi_u(i,m+1) = pi/4 + temp*pi/2;
    end
    
    cvx_begin quiet
    variable w
    variable X(N,2^B) binary
    minimize w
    subject to
    w >= abs(imag(H*X*v.*exp(-1i*phi_u(:,m+1))))-real(H*X*v.*exp(-1i*phi_u(:,m+1)));
    for i = 1:N
        sum(X(i,:)) == 1;
    end
    cvx_end
    X_o(:,m+1) = X*v;
end
t_o = min(min((real(H*X_o.*exp(-1i.*phi_u))-abs(imag(H*X_o.*exp(-1i.*phi_u))))/sqrt(2)));