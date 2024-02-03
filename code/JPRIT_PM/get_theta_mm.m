% Solve for the problem (43) to obtain theta.
% This is used in the paper: R. Liu, M. Li, Q. Liu, A. L. Swindlehurst, and Q. Wu,“Intelligent reflecting surface based passive information transmission: A symbol-level precoding approach,” IEEE Trans. Veh. Technol., vol. 70, no. 7, pp. 6735-6749, Jul. 2021.
% Download this paper at: https://ieeexplore.ieee.org/document/9435988
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-02

function [theta] = get_theta_mm(a,B,K,theta)

[nrow,~] = size(a);
[~,N] = size(B);

Bi = zeros(nrow,N);
Ci = zeros(nrow,N);
Ui = zeros(nrow,N);
Vi = zeros(nrow,N);
wi = zeros(nrow,1);
zi = zeros(nrow,1);


phi = pi/4;
for iar = 1:1:2*K*4^K
    Bi(iar,:) = cos(phi).*imag(B(iar,:)) - real(B(iar,:)).*sin(phi);
    Ci(iar,:) = cos(phi).*real(B(iar,:)) + imag(B(iar,:)).*sin(phi);
    Ui(iar,:) = -cos(phi).*imag(B(iar,:)) - real(B(iar,:)).*sin(phi);
    Vi(iar,:) = imag(B(iar,:)).*sin(phi) - cos(phi).*real(B(iar,:));
    wi(iar) = imag(a(iar))*cos(phi) - real(a(iar))*sin(phi);
    zi(iar) = -imag(a(iar))*cos(phi) - real(a(iar))*sin(phi);
end
phi = pi/2;
for iar = 2*K*4^K+1:nrow
    Bi(iar,:) = cos(phi).*imag(B(iar,:)) - real(B(iar,:)).*sin(phi);
    Ci(iar,:) = cos(phi).*real(B(iar,:)) + imag(B(iar,:)).*sin(phi);
    Ui(iar,:) = -cos(phi).*imag(B(iar,:)) - real(B(iar,:)).*sin(phi);
    Vi(iar,:) = imag(B(iar,:)).*sin(phi) - cos(phi).*real(B(iar,:));
    wi(iar) = imag(a(iar))*cos(phi) - real(a(iar))*sin(phi);
    zi(iar) = -imag(a(iar))*cos(phi) - real(a(iar))*sin(phi);
end

% Create the problem structure.
manifold = obliquefactory(2, N);
problem.M = manifold;

% Define the problem cost function and its gradient.
problem.cost = @cost;
    function h = cost(x)
        f = Bi*(x(1,:).') + Ci*(x(2,:).') + wi;
        g = Ui*(x(1,:).') + Vi*(x(2,:).') + zi;
        h = epsl*log( sum( exp(f/epsl)+exp(g/epsl) ) );
    end

problem.grad = @(x) problem.M.egrad2rgrad(x,egrad(x));
    function eh = egrad(x)
        f = Bi*(x(1,:).') + Ci*(x(2,:).') + wi;
        g = Ui*(x(1,:).') + Vi*(x(2,:).') + zi;
        h1 = sum(repmat(exp(f/epsl),1,N).*Bi) + sum(repmat(exp(g/epsl),1,N).*Ui);
        h2 = sum(repmat(exp(f/epsl),1,N).*Ci) + sum(repmat(exp(g/epsl),1,N).*Vi);
        h3 = sum( exp(f/epsl)+exp(g/epsl) );
        eh = [h1;h2]./h3;
    end

% Execute the optimization
options.tolgradnorm = 1e-10;
options.maxiter = 500;
options.minstepsize = 1e-20;
options.verbosity = 0;
epsl = 0.005;
tag = 0;
while 1
    [x,~,info] = conjugategradient(problem,[real(theta.');imag(theta.')],options);

    if isnan(x)
        if tag == 3
            tag = 0;
            epsl = epsl + 0.001;
        else
            tag = tag + 1;
            theta = exp(1i*2*pi*rand(N,1));
        end
    else
        theta = (x(1,:)+1i*x(2,:)).';
        break;
    end
end

end

