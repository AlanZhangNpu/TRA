function result = MCS( problem, par )
% Monto Carlo Simulation 
% for time-variant reliability analysis

% setting parameters
if nargin == 3
    n_MCS = 1e6;
    n_time_instant = 400;
else
    n_MCS = par.n_MCS;
    n_time_instant = par.n_time_instant;
end
variable_table = problem.variable_table;
performanceFunc = problem.performanceFunc;

% identify the number of static and dynamic RVs
n_RV = size(variable_table,1);
n_static_RV = n_RV;
for rv_id = 1:n_RV
    if ~isempty(variable_table{rv_id,end})
        n_static_RV = rv_id - 1;
        break;
    end
end
n_dynamic_RV = n_RV - n_static_RV;
time_series = linspace(problem.Ts, problem.Te, n_time_instant);

% run simulation
first_failure_list = inf*ones(n_MCS,1);
batch_size = 5e5;
batch_num = floor(n_MCS/batch_size);
for batch_id = 1:batch_num
    start_id = (batch_id-1) * batch_size + 1;
    stop_id = batch_id * batch_size;
    first_failure_list(start_id:stop_id) = batch(batch_size);
    
    disp([datestr(clock) '    ' num2str(stop_id/n_MCS*100) '%']);
end
left_num = n_MCS - batch_size*batch_num;
if left_num > 0
    start_id = batch_num * batch_size + 1;
    first_failure_list(start_id:end) = batch(left_num);    
    disp([datestr(clock) '    100%']);
end

% calculate the evolution of pf
pf_history = zeros(1,n_time_instant);
for time_point_id = 1:n_time_instant
    failure_num = sum(first_failure_list <= time_point_id);
    pf_history(time_point_id) = failure_num/n_MCS;
end

result.pf = pf_history(end);
result.pf_history = [time_series; pf_history];
result.funccount = n_MCS*n_time_instant;
return;

    function first_failure_list = batch(n_MCS)
        % realization of the static RVs
        static_RVs = zeros(n_MCS, n_static_RV);
        for rv_id = 1:n_static_RV
            static_RVs(:,rv_id) = ...
                GenerateRV( ...
                variable_table{rv_id,1}, ...
                variable_table{rv_id,2}, ...
                variable_table{rv_id,3}, ...
                n_MCS);
        end
        
        % realization of the dynamic RVs
        dynamic_RVs = cell(1, n_dynamic_RV);
        for dynamic_rv_id = 1 : n_dynamic_RV
            dynamic_RVs{dynamic_rv_id} = ...
                GenerateGP( ...
                variable_table{n_static_RV + dynamic_rv_id,1}, ...
                variable_table{n_static_RV + dynamic_rv_id,2}, ...
                variable_table{n_static_RV + dynamic_rv_id,3}, ...
                variable_table{n_static_RV + dynamic_rv_id,4}, ...
                time_series, n_MCS);
        end
        
        % run simulation
        first_failure_list = inf*ones(n_MCS,1);
        for simulation_id = 1:n_MCS            
            static_part = static_RVs(simulation_id, :);
            dynamic_part = zeros(n_dynamic_RV, n_time_instant);
            for dynamic_rv_id = 1 : n_dynamic_RV
                dynamic_part(dynamic_rv_id,:) = dynamic_RVs{dynamic_rv_id}(simulation_id, :);
            end
            
            G = performanceFunc(time_series, static_part, dynamic_part);
            first_failure_list(simulation_id) = calFirstFailure(G);
        end        
    end

    function first_failure = calFirstFailure(G)
        label = G >= 0;
        first_failure = find(label,1);
        if isempty(first_failure)
            first_failure = inf;
        end
    end

end

