function Pf = CalculatePf(MPPList, time_series, ACF, inverse)
% calculate the probability of failure

% step 1: calculate ALPHA and BETA
[n_MPP, n_var] = size(MPPList);
BETA = zeros(1,n_MPP);
ALPHA = zeros(size(MPPList));
for i=1:n_MPP
    BETA(i) = norm(MPPList(i,:));
    ALPHA(i,:) = MPPList(i,:) / BETA(i);
end

% step 2: calculate the covariance matrix SIGMA
SIGMA = ones(n_MPP,n_MPP);
for i=1:n_MPP
    for j=1:n_MPP
        if i ~= j
            t1 = time_series(i);
            t2 = time_series(j);
            alpha1 = ALPHA(i,:);
            alpha2 = ALPHA(j,:);
            SIGMA(i,j) = calResponseCC(t1, alpha1, t2, alpha2);
        end
    end
end

% step 3:
if nargin == 4 && inverse
    BETA = -BETA;
end
Pf = 1 - mymvncdf(BETA, SIGMA);
return;

    function rho = calResponseCC(t1, alpha1, t2, alpha2)
        C = eye(n_var); % identify matrix
        for ii = 1:n_var
            if ~isempty(ACF{ii})
                rho_x = ACF{ii}(t1, t2);
                rho_u = rho_x;              % NOTE: only for Gaussin process
                C(ii, ii) = rho_u;
            end
        end
        rho = alpha1 * C * alpha2';
    end
end


%% CDF for multivariant normal distribution through MCS
function p = mymvncdf(beta, SIGMA)
% the following command is not reliable when 
% the dimension is high
% p = mvncdf(beta, zeros(1,length(beta)), SIGMA);

p_list = zeros(10,1);
for trial_id = 1:length(p_list)
    p_list(trial_id) = MCS(beta, SIGMA);
end
p = mean(p_list);
return;

    function p = MCS(beta, SIGMA)
        n_total = 1e6;
        sample = mvnrnd(zeros(1,length(beta)), SIGMA, n_total);
        label = sum(sample < beta,2) ~= length(beta);
        n_fail = sum(label);
        n_rare_event = min(n_fail, n_total - n_fail);
        while n_rare_event < 100 && n_total < 1e7
            if n_rare_event == 0
                n_rare_event = 1;
            end
            n_additional_MCS = ceil(n_total * (150/n_rare_event - 1));
            
            % considering the memory limit
            n_additional_MCS = min(n_additional_MCS, 1e6);
            
            sample = mvnrnd(zeros(1,length(beta)), SIGMA, n_additional_MCS);
            label = sum(sample < beta,2) ~= length(beta);
            
            n_fail = n_fail + sum(label);
            n_total = n_total + n_additional_MCS;
            n_rare_event = min(n_fail, n_total - n_fail);
        end
        p = (n_total - n_fail)/n_total;
    end
end