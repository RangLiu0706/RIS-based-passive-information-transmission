% Solve for the problem (7).
% This is used in the paper: R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02
% Inputs: H: the channel;
% Outputs: X_CE: RIS reflecting coefficients;
%          t: the objective value
function [X_CE,t] = get_X_CE(H)
[K,N] = size(H);
phi_u = zeros(K,4^K);
X_CE = zeros(N,4^K);
A = zeros(K,N);
for m = 0:1:4^K-1    
    s = dec2bin(m,2*K);
    for i = 1:1:K
        temp = 2*(s(2*i-1)-48)+s(2*i)-48;
        phi_u(i,m+1) = pi/4 + temp*pi/2;
        A(i,:) = H(i,:).*exp(-1i*phi_u(i,m+1));
    end    
    
    X_CE(:,m+1) = oblique_CE(A.*10^5);         
end
t = min(min((real(H*X_CE.*exp(-1i.*phi_u))-abs(imag(H*X_CE.*exp(-1i.*phi_u))))/sqrt(2)));