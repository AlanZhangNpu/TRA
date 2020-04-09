function result = AMPPT( problem, option )
% AMPPT
% Efficient Time-variant Reliability Analysis
% through Approximating the MPP trajectory
% 
% Structure of the arguments
% problem:
%   problem.variable_table;
%   problem.Ts;
%   problem.Te;
%   problem.performanceFunc;
% 
% option:
%   option.approximation_initial_point_num
%   option.approximation_tolerance
%   option.discretization_node_num
% 
% auther: zyw (ywzhang@nwpu.edu.cn)


if ~isfield(option, 'inverse')
    option.inverse = 0;
end
if ~isfield(option, 'approximation_initial_point_num')
    option.approximation_initial_point_num = 3;
end
if ~isfield(option, 'approximation_tolerance')
    option.approximation_tolerance = 1e-4;
end

[MPPT_model, statistics] = createMPTModel();
pf_history = calReliability(MPPT_model);
for k = 1:size(pf_history,2)
    pf_history(2,k) = max(pf_history(2,1:k));
end

result.pf = pf_history(2,end);
result.pf_history = pf_history;
result.funccount = sum(statistics(:,2));
result.statistics = statistics;
return;

    function model = train(time_series, mpp_list)        
        [mpp_list_n,settings] = mapminmax(mpp_list');
        mpp_list_n = mpp_list_n';
        
        n_y_dim = size(mpp_list,2);
        total_model = cell(1,n_y_dim);
        theta = 1; lob = 0.2; upb = 2;
        for ii = 1:n_y_dim
            total_model{ii} = dacefit( ...
                time_series, mpp_list_n(:,ii), @regpoly0, @corrgauss, theta, lob, upb);
        end
        model = @predict;
        
        function [y,mse] = predict(x)
            for jj = 1:n_y_dim
                [y(:,jj), mse(:,jj)] = predictor(x, total_model{jj});
            end
            y = mapminmax.reverse(y',settings)';
        end
    end

    function [MPPT_model, statistics] = createMPTModel()
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
        
        fprintf('%8s%16s%16s%16s%16s\n', ...
            'mpp_num', 'funccount', 'total funccount', 't_next', 'max_mse');
        
        % step 1: initial sampling: generate equally spaced samples
        time_series = linspace(problem.Ts, problem.Te, option.approximation_initial_point_num);
        mpp_list = zeros(option.approximation_initial_point_num, size(variable_table,1));
        funccout_WS = zeros(option.approximation_initial_point_num,1);
        funccount_NWS = zeros(option.approximation_initial_point_num,1);
        error = zeros(option.approximation_initial_point_num,1)*NaN;
        for id = 1:option.approximation_initial_point_num
            start_point_u = u0;
            if id ~= 1
                start_point_u = mpp_list(id-1,:);
            end
            
            variable_table_copy = copy_table(variable_table, time_series(id));
            
            [mpp_list(id,:), funccout_WS(id)] = SearchMPP( ...
                @(x) performanceFunc(time_series(id), x(1:n_static_RV), x(n_static_RV+1:end)'), ...
                variable_table_copy, ...
                start_point_u);
            
            if id < option.approximation_initial_point_num
                fprintf('%8g%16g%16g\n', id, funccout_WS(id), sum(funccout_WS));
            end
        end
        
        % step 2: adaptive sampling to approximate the Most Probable Curve (MPC)
        while 1
            % update the Kriging model
            MPPT_model = train(time_series', mpp_list);
                        
            % adaptive sampling
            subproblem.solver = @(x) solver(x, MPPT_model);
            subproblem.bound = [problem.Ts, problem.Te];
            subproblem.popSize = 100;
            subproblem.maxEvaluation = 1e8;
            [optimum, ~] = DE_Canonical(subproblem);
            t_next = optimum.best_point;
            max_mse = -optimum.best_fitness;            
            error(id) = max_mse;
            
            fprintf('%8g%16g%16g%16g%16g\n', ...
                size(mpp_list,1), funccout_WS(id), sum(funccout_WS), t_next, max_mse);
            
            if ~all(time_series - t_next)
                break;
            end
            
            if max_mse <= option.approximation_tolerance
                break;
            end
            
            % run MPP search at the new generated time point
            id = id + 1;
            time_series(id) = t_next;            
            variable_table_copy = copy_table(variable_table, t_next);
            
            % warm-start: predict the MPP at t_next
            start_point_u = MPPT_model(t_next);
            
            [mpp_list(id,:), funccout_WS(id)] = SearchMPP( ...
                @(x) performanceFunc(t_next, x(1:n_static_RV), x(n_static_RV+1:end)'), ...
                variable_table_copy, ...
                start_point_u);
        end
        
%         statistics = [time_series' funccout_WS funccount_NWS mpp_list error];
        statistics = [time_series' funccout_WS mpp_list error];
        
        function result = solver(x, MPPT_model)
            [~, mse] = MPPT_model(x);
            result = zeros(size(x,1),1);
            for individual_id = 1:size(x,1)
                result(individual_id) = -mean(mse(individual_id,:));
            end
        end        
    end

    function pf_history = calReliability(MPPT_model)
        % step 3: discretization of the approximated MPPT
        time_series = linspace(problem.Ts, problem.Te, option.discretization_node_num);
        mpp_list = MPPT_model(time_series');
        
        % step 4: calculate the probability of failure through 
        % spectral decomposition and MCS
        n_var = size(problem.variable_table,1);
        ACF = cell(n_var, 1);
        for var_id = 1:n_var
            ACF{var_id} =  problem.variable_table{var_id, 4};
        end
        
        fprintf('%8s%16s%16s\n', 'id', 'time', 'pf');
        
        n_time_point_export = min(40, option.discretization_node_num);
        time_points = round(linspace(1, option.discretization_node_num, n_time_point_export));
        pf_history = zeros(2, n_time_point_export);
        for jj = 1:n_time_point_export
            time_point_id = time_points(jj);
            pf_history(1, jj) = time_series(time_point_id);
            pf_history(2, jj) = CalculatePf(mpp_list(1:time_point_id,:), ...
                time_series(1:time_point_id),ACF,option.inverse);
            
            fprintf('%8d%16f%16f\n', jj, pf_history(1, jj), pf_history(2, jj));
        end
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
