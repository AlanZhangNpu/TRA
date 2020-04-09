function [MPP, funccout] = SearchMPP(performanceFunc, variable_table, x0)
% search MPP

ub = 10 * ones(1,size(variable_table,1));
options = optimoptions('fmincon','Display','none');
[MPP,fval,exitflag,output] = fmincon(@obj,x0,[],[],[],[],-ub,ub,@con,options);
if exitflag ~= 1 && output.constrviolation > 1e-5 && fval< -norminv(1e-4)
    warning(['Failed to search the MPP: ' output.message]);
end
funccout = output.funcCount;

    function f = obj(u)
        f = norm(u);
    end

    function [c,ceq] = con(u)
        c = [];
        x = NatafTransformation(u, variable_table, -1);
        ceq = performanceFunc(x);
    end
end

