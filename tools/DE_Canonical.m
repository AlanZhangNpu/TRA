function [best, result] = DE_Canonical(problem)
% the canonical DE algorithm
% the parameter "problem" must have the following members:
%     problem.solver.obj
%     problem.solver.con
%     problem.bound
%     problem.popSize
%     problem.maxEvaluation

F = 0.9;
Cr = 0.4;
MAX_STALLED_ITERATIONS = 40;

x_dim = size(problem.bound,1);

pop = [];
fitness = [];
constraint = [];

best.best_point = NaN;
best.best_fitness = NaN;
best.best_constraint = NaN;

funccount = 0;
stall_iterations = 0;
is_continue = 1;

% tic();
while is_continue == 1
    globalSearch();
    if funccount >= problem.maxEvaluation
        is_continue = -1;
    end
    if stall_iterations >= MAX_STALLED_ITERATIONS
        is_continue = -2;
    end
end
% result.elapsed_time = toc();
result.funccount = funccount;
if is_continue == -1
    result.msg = ['The user-specified maxEvaluation(' num2str(problem.maxEvaluation) ') is reached.'];
elseif is_continue == -2
    result.msg = ['MAX_STALLED_ITERATIONS(' num2str(MAX_STALLED_ITERATIONS) ') is reached.'];
end


%% subfunctions for global search
    function globalSearch()
        if isempty(pop)
            % initialization
            pop = DOE_latin( problem.bound, problem.popSize );
            res = problem.solver(pop);
            fitness = res(:,1);
            constraint = res(:,2:end);
        else
            % obtain the next population
            newpop = move(pop, problem.bound);
            
            res = problem.solver(newpop);
            newfitness = res(:,1);
            newContraint = res(:,2:end);
            
            for i = 1:problem.popSize
                if isempty(constraint)
                    superior = comparePoints(fitness(i), [], newfitness(i), []);
                else
                    superior = comparePoints(fitness(i), constraint(i,:), newfitness(i), newContraint(i,:));
                end
                if superior == 1
                    pop(i,:) = newpop(i,:);
                    fitness(i) = newfitness(i);
                    if ~isempty(constraint)
                        constraint(i,:) = newContraint(i,:);
                    end
                end
            end            
        end
        
        funccount = funccount + problem.popSize;
        
        % archive
        best_id = findBestIndividual(fitness, constraint); % findBestFeasibleIndividual
        if isnan(best.best_fitness) || comparePoints(best.best_fitness, best.best_constraint, ...
                fitness(best_id,:), constraint(best_id,:)) == 1
            best.best_point = pop(best_id,:);
            best.best_fitness = fitness(best_id,:);
            best.best_constraint = constraint(best_id,:);
            stall_iterations = 0;
        else
            stall_iterations = stall_iterations + 1;
        end
    end

    function U = move(X, bound)
        F = 0.8; Cr = 0.4;
        popSize = problem.popSize;
        V = X;
        for ii = 1:popSize
            squence = 1:popSize;
            squence(ii) = [];
            r = randperm(popSize-1,3);
            r1 = squence(r(1));
            r2 = squence(r(2));
            r3 = squence(r(3));
            V(ii, :) = X(r1,:) + F * (X(r2,:) - X(r3,:));
        end
        indexCross = logical(rand(popSize, x_dim) <= Cr | repmat(1 : x_dim, popSize, 1) == repmat(randi(x_dim, [popSize, 1]), 1, x_dim));
        U = V .* indexCross + X .* (1 - indexCross);
        U = box(U, bound);
        
        function pop = box(pop, bound)
            for popIndex = 1 : popSize
                for dimIndex = 1:x_dim
                    if pop(popIndex,dimIndex) < bound(dimIndex, 1)
                        pop(popIndex,dimIndex) = bound(dimIndex, 1);
                    end
                    if pop(popIndex,dimIndex) > bound(dimIndex, 2)
                        pop(popIndex,dimIndex) = bound(dimIndex, 2);
                    end
                end
            end
        end
    end

    function best_id = findBestIndividual(fitness, constraint)
        best_id = 1;
        for i = 2 : problem.popSize
            superior = comparePoints(fitness(best_id,:), constraint(best_id,:),fitness(i,:), constraint(i,:));
            if superior == 1
                best_id = i;
            end
        end
    end

    function is_feasible = isFeasible(constraint)
        constraint_num = size(constraint,2);
        if constraint_num == 0 || sum(constraint <= 0) == constraint_num
            is_feasible = 1;
        else
            is_feasible = 0;
        end
    end

    function [ x ] = DOE_latin( bound, n )
        % Latin Hypercube Sampling 
        % for generating the initial population
        
        dim = size(bound, 1); % dimension
        x = lhsdesign(n, dim, 'criterion','maximin');
        for i = 1 : dim
            down = bound(i, 1);
            up = bound(i, 2);
            x(:, i) = x(:, i) .* (up - down) + down;
        end
    end

    function moveon = comparePoints(fitness1, contraint1, fitness2, contraint2)
        % superior == 0: the first point is better
        % superior == 1: the second point is better
        
        v1 = calV(contraint1);
        v2 = calV(contraint2);
        
        if v1 == 0 && v2 == 0
            moveon = compareFitness(fitness1, fitness2);
            return;
        end
        moveon = compareFitness(v1, v2);
        return;
        
        function v = calV(contraint)
            contraint(1,contraint<=0) = 0;
            v = sum(contraint);
        end
        
        function moveon = compareFitness(f1, f2)
            if f1 > f2
                moveon = 1;
            else
                moveon = 0;
            end
        end
        
    end
end