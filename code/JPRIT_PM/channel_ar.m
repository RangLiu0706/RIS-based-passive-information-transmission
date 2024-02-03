% Generate the channel.
% This is used in the paper: R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02
function [H]=channel_ar(M,N)

ind_Tx = (0:1:M-1)';
ind_Rx = (0:1:N-1)';
AoD = pi*rand - pi/2;
AoA = pi*rand - pi/2;

Amh = sqrt(1/N) * exp( 1j*pi*ind_Rx*sin(AoA) );
Abh = sqrt(1/M) * exp( 1j*pi*ind_Tx*sin(AoD) );

H = sqrt(M*N)*Amh*Abh';
end
