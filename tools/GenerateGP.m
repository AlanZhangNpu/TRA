function [ realization ] = GenerateGP( type, meanFunc, stdFunc, autocorrelationCoefFunc, time_series, n )
% generate Gaussian Processes

is_stationary = strcmp(type, 'stationary');
nnode = length(time_series);

% calculate the mean value vector
if is_stationary
    mu = meanFunc * ones(size(time_series));
else
    mu = zeros(size(time_series));
    for ii=1:nnode
        mu(ii) = meanFunc(time_series(ii));
    end
end

% calculate the covariance matrix
sigma = zeros(nnode, nnode);
for ii = 1 : nnode
    t1 = time_series(ii);
    if is_stationary
        sigma(ii,:) = stdFunc .* stdFunc .* autocorrelationCoefFunc(t1, time_series);
    else
        sigma(ii,:) = stdFunc(t1) .* stdFunc(time_series) .* autocorrelationCoefFunc(t1, time_series);
    end
end

% generate the normally distributed variables
realization = mvnrnd(mu, sigma, n);

end