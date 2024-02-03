% Calculated the SER.
% This is used in the paper: R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02
% Inputs: R: the received signals; S: the tranmistted symbols; 
%         omega: the modulation order;
% Outputs: ser: the SER.
function ser = get_SER(R,S,omega)

[K,N_s] = size(R);
ser = zeros(1,K);
for ik = 1:1:K
    mpsk = exp( 1i*(pi/omega(ik)+2*pi/omega(ik)*(0:1:omega(ik)-1)) );
    Distance = abs(repmat(R(ik,:),omega(ik),1)-repmat(mpsk.',1,N_s));
    [~,R(ik,:)] = min(Distance); 
    [~,ser(ik)] = symerr(R(ik,:),S(ik,:));
end




