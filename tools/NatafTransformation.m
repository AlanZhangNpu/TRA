function u = NatafTransformation( x, variable_table, ~ )
%Nataf transformation
% transform the input variables into 
% the standard normal space (u space)

if nargin == 2
    forward = 1;
else
    forward = 0;
end
u = zeros(size(x));
for variable_id = 1:size(variable_table,1)
    if forward
        func = @Nataf;
    else
        func = @invNataf;
    end
    u(:, variable_id) = func(x(:, variable_id), ...
        variable_table{variable_id, 1}, ...
        variable_table{variable_id, 2}, ...
        variable_table{variable_id, 3});
end
return;
end

function u = Nataf(x, distribution, mu, sigma)
mid = mycdf(x, distribution, mu, sigma);
u = mycdf(mid, 'normalinv');
end

function x = invNataf(u, distribution, mu, sigma)
mid = mycdf(u);
x = mycdf(mid, [distribution 'inv'], mu, sigma);
end


%% CDF
function p = mycdf(x, distribution, mu, sigma)
if nargin == 1
    distribution = 'normal';
    mu = 0;
    sigma = 1;
end
if nargin == 2
    mu = 0;
    sigma = 1;
end

if strcmp(distribution, 'normal') || strcmp(distribution, 'stationary')
    p = normcdf(x, mu, sigma);
elseif strcmp(distribution, 'normalinv') || strcmp(distribution, 'stationaryinv')
    p = norminv(x, mu, sigma);
    
elseif strcmp(distribution, 'lognormal')
    [mu, sigma] = lognormtrans(mu, sigma);
    p = logncdf(x, mu, sigma);
elseif strcmp(distribution, 'lognormalinv')
    [mu, sigma] = lognormtrans(mu, sigma);
    p = logninv(x, mu, sigma);
    
elseif strcmp(distribution, 'gumbel')
    [u, alpha] = gumbeltrans(mu, sigma);
    p = evcdf(x, u, alpha);
elseif strcmp(distribution, 'gumbelinv')
    [u, alpha] = gumbeltrans(mu, sigma);
    p = evinv(x, u, alpha);
    
else
    error(['Unkonwn distribution: ' distribution]);
end
end

function [u, alpha] = gumbeltrans(mean_value, std_value)
gama = -psi(1);
alpha = sqrt(6)*std_value/pi;
u = gama*alpha + mean_value;
end

function [MU, SIGMA] = lognormtrans(mean_value, std_value)
% transformation for lognormal distribution
V = std_value^2;       % variance
MU = log(mean_value^2 / sqrt(V+mean_value^2));
SIGMA = sqrt(log(V/mean_value^2 + 1));
end
