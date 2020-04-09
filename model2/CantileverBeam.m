function g = CantileverBeam( t, random_variables, stochastic_processes )
% performance function of the cantilever beam

d   = random_variables(:,1);
h   = random_variables(:,2);
R0  = random_variables(:,3);
F2  = random_variables(:,4);
P   = random_variables(:,5);

F1   = stochastic_processes(1,:);
T_t  = stochastic_processes(2,:);

L1 = 0.06;
L2 = 0.12;
theta1 = deg2rad(10);
theta2 = deg2rad(5);

A = pi/4  .*(d.^2 - (d-2*h).^2);
I = pi/64 .*(d.^4 - (d-2*h).^4);
M = F1.*cos(theta1).*L1 + F2.*cos(theta2).*L2;
sigma_t = (F1.*sin(theta1) + F2.*sin(theta2) + P) ./ A + M.*d ./ (2*I);
tao_t = T_t .* d ./ (4*I);
sigma_max_t = sqrt(sigma_t.^2 + 3*tao_t.^2);
R_t = R0.*(1 - 0.01*t);
g = sigma_max_t - R_t;

    function rad = deg2rad(deg)
        rad = deg ./ 180 .* pi;
    end

end