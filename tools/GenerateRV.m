function res = GenerateRV( method, mean_value, std_value, row_count, col_count )
% generate random variables

if nargin == 3
    row_count = 1e6;
    col_count = 1;
end

if nargin == 4
    col_count = 1;
end

if strcmp(method, 'normal')
    res = normrnd(mean_value,std_value,row_count,col_count);
    
elseif strcmp(method, 'lognormal')
    V = std_value^2;
    MU = log(mean_value^2 / sqrt(V+mean_value^2));
    SIGMA = sqrt(log(V/mean_value^2 + 1));
    res = lognrnd(MU, SIGMA, row_count,col_count);
    
elseif strcmp(method, 'gumbel')        
    gama = -psi(1);
    alpha = sqrt(6)*std_value/pi;
    u = gama*alpha + mean_value;        
    res = evrnd(u,alpha,row_count, col_count);
end

end

