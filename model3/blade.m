function g = blade(t, random_variable, stochastic_process)
l1 = random_variable(:,1);
t1 = random_variable(:,2);
t2 = random_variable(:,3);
e_allowable = random_variable(:,4);

v_t = stochastic_process;

rho = 1e3;
Cm = 0.3422;

M_flap = 0.5 * rho .* v_t.^2 .* Cm;
E = 14e9;
I = 2/3*l1.*(t1.^3 - t2.^3);

a = M_flap.*t1 ./ (E.*I);
g = e_allowable - a;
g = -g;
end


