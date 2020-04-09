function result = TDTRA( problem, option )
% Time-discritization based time-variant
% reliability analysis (TDTRA)
% 
% Structure of the arguments
% problem:
%   problem.variable_table;
%   problem.Ts;
%   problem.Te;
%   problem.performanceFunc;
% 
% option:
%   option.inverse
%   option.time_node_num
% 
% auther: zyw (ywzhang@nwpu.edu.cn)

if ~isfield(option, 'inverse')
    option.inverse = 0;
end

% pf_history
[pf_history, statistics] = analyze(problem);

for ii=1:size(pf_history,2)
    pf_history(2,ii) = max(pf_history(2,1:ii));
end

result.pf = pf_history(2,end);
result.funccount = sum(statistics(:,2));
result.pf_history = pf_history;
result.statistics = statistics;


    function [pf_history, statistics] = analyze(problem)
        
        % parameters
        variable_table = problem.variable_table;
        performanceFunc = problem.performanceFunc;
        n_var = size(variable_table,1);
        n_static_RV = n_var;
        for rv_id = 1:n_var
            if ~isempty(variable_table{rv_id,end})
                n_static_RV = rv_id - 1;
                break;
            end
        end
        
        % start point of MPP search in the original space
        u0 = zeros(1, n_var);
        
        fprintf('%4s%16s%16s%16s%16s\n', ...
            'id', 'time', 'RI', 'funccount', 'total funccount');
        
        % step 1: time discretization and sampling
        time_series = linspace(problem.Ts, problem.Te, option.time_node_num);
        mpp_list = zeros(option.time_node_num, size(variable_table,1));
        funccout_list = zeros(option.time_node_num,1);
        for id = 1:option.time_node_num
            variable_table_copy = copy_table(variable_table, time_series(id));
            [mpp_list(id,:), funccout_list(id)] = SearchMPP( ...
                @(x) performanceFunc(time_series(id), x(1:n_static_RV), x(n_static_RV+1:end)'), ...
                variable_table_copy, ...
                u0);
            fprintf('%4g%16g%16g%16g%16g\n', ...
                id, time_series(id), norm(mpp_list(id,:)), funccout_list(id), sum(funccout_list));
        end
        
        % step 2: calculate the probability of failure through series expansion of
        % the stochastic process
        n_var = size(problem.variable_table,1);
        ACF = cell(n_var, 1);
        for var_id = 1:n_var
            ACF{var_id} =  problem.variable_table{var_id, 4};
        end
        
        fprintf('%8s%16s%16s\n', 'id', 'time', 'pf');
        
        n_time_point_export = min(40, option.time_node_num); % 10
        pf_history = zeros(2, n_time_point_export);
        time_points = round(linspace(1, option.time_node_num, n_time_point_export));
        for jj = 1:n_time_point_export
            time_point_id = time_points(jj);
            
            pf_history(1, jj) = time_series(time_point_id);
            pf_history(2, jj) = CalculatePf( ...
                mpp_list(1:time_point_id,:), ...
                time_series(1:time_point_id), ...
                ACF, ...
                option.inverse);
            fprintf('%8d%16f%16f\n', jj, pf_history(1, jj), pf_history(2, jj));
        end
        
        statistics = [time_series' funccout_list mpp_list];
    end

    function variable_table_copy = copy_table(vt, t)
        variable_table_copy = vt;
        for rv_id = 1:size(vt,1)
            if isa(variable_table_copy{rv_id,2},'function_handle')
                variable_table_copy{rv_id,2} = variable_table_copy{rv_id,2}(t);
            end
            if isa(variable_table_copy{rv_id,3},'function_handle')
                variable_table_copy{rv_id,3} = variable_table_copy{rv_id,3}(t);
            end
            if ~isempty(variable_table_copy{rv_id,4})
                variable_table_copy{rv_id,1} = 'normal';
            end
        end
    end
end