% Solve for the problem (11).
% This is used in the paper: R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02
function [x_ce] = oblique_CE(A)

[K,N] = size(A);
phi = pi/4; % QPSK
B = zeros(K,N);
C = zeros(K,N);
U = zeros(K,N);
V = zeros(K,N);
for i = 1:1:K
    B(i,:) = cos(phi).*imag(A(i,:)) - real(A(i,:)).*sin(phi);
    C(i,:) = cos(phi).*real(A(i,:)) + imag(A(i,:)).*sin(phi);
    U(i,:) = -cos(phi).*imag(A(i,:)) - real(A(i,:)).*sin(phi);
    V(i,:) = imag(A(i,:)).*sin(phi) - cos(phi).*real(A(i,:));
end

% Create the problem structure.
manifold = obliquefactory(2, N);
problem.M = manifold;
% Define the problem cost function and its gradient.

problem.cost = @cost;
    function f = cost(x)
        f1 = B*(x(1,:).') + C*(x(2,:).');
        f2 = U*(x(1,:).') + V*(x(2,:).');
        f_temp = exp(f1./epsl) + exp(f2./epsl);
        f = epsl*log(sum(f_temp));
    end
problem.grad = @(x) problem.M.egrad2rgrad(x,egrad(x));
    function g = egrad(x)
        f1 = B*(x(1,:).') + C*(x(2,:).');
        f2 = U*(x(1,:).') + V*(x(2,:).');
        f3 = sum( exp(f1./epsl) + exp(f2./epsl) );
        
        f4 = sum(B.*repmat( exp(f1./epsl),1,N ) ) + sum(U.*repmat( exp(f2./epsl),1,N ) );
        f5 = sum(C.*repmat( exp(f1./epsl),1,N ) ) + sum(V.*repmat( exp(f2./epsl),1,N ) );

        g = [f4;f5]./f3;
    end

% Execute the optimization
options.tolgradnorm = 1e-6;
options.maxiter = 500;
options.verbosity = 0;

epsl = 0.1;
[x,~,info,~] = conjugategradient(problem,[],options);
x_ce = (x(1,:)+1i*x(2,:)).';
end