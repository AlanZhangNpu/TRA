function [ G ] = SteelBeam( t, random_variable, stochastic_processes )
%STEELBEAM 

% b0, h0, sigmay, Ft
b0 = random_variable(:,1);
h0 = random_variable(:,2);
sigmay = random_variable(:,3);
Ft = stochastic_processes;

k=5e-5;
L=5;
rho=78.5e3;

bt = b0 - 2*k*t;
ht = h0 - 2*k*t;

term1 = bt .* ht .^ 2 .* sigmay ./ 4;
term2 = Ft.*L./4 + rho .* b0 .* h0 .* L.^2 ./ 8;
G = term2 - term1;
end