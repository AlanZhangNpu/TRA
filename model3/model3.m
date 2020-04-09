function model3()

rng('shuffle');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% problem definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
problem.performanceFunc = @blade;
problem.Ts = 0;
problem.Te = 12;
problem.variable_table = {
    % distribution      mean     	std                 autocorrelationCoef
    % l1, t1, t2, e_allowable
    'normal',           0.22,       0.0022,         [];
    'normal',           0.025,      0.00025,        [];
    'normal',           0.019,      0.00019,     	[];
    'normal',           0.025,   	0.00025,      	[];
    % v(t)
    'nonstationary',    @meanFunc, 	@stdFunc,    	@autoCoefFunc;
    };

    function m = meanFunc(t)
        a_m = [3.815,   2.528,  1.176,  -0.07856];
        b_m = [0.2895,  0.5887, 0.7619, 2.183];
        c_m = [-0.2668, 0.9651, 3.116,  -3.161];
        temp = a_m .* sin(b_m .* t + c_m);
        m = sum(temp);
    end

    function s = stdFunc(t)
        t = t';
        a_s = [0.7382,  1.013,  1.875,  1.283];
        b_s = [6.456,   4.075,  9.913,  1.035];
        c_s = [0.9193,  1.561,  6.959,  2.237];
        temp = a_s .* exp(-((t-b_s)./c_s).^2);
        s = sum(temp,2)';
    end
    function rho = autoCoefFunc(t1, t2)
        rho = cos(2*pi*(t2-t1));
    end

solve(problem);
analyzeResults();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function solve(problem)
% MCS
MCS_par.n_MCS = 1e6;
MCS_par.n_time_instant = 300;
result.MCS = MCS(problem, MCS_par);

% ADMPT
option.discretization_node_num = 200;
result.AMPPT = AMPPT( problem, option );

% TDTRA
option.time_node_num = 100;  result.TDTRA100 = TDTRA( problem, option );
option.time_node_num = 150;  result.TDTRA150 = TDTRA( problem, option );
option.time_node_num = 200;  result.TDTRA200 = TDTRA( problem, option );
save('result.mat', 'result');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% analyze the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analyzeResults()
result = [];
load('result.mat');

accuracy = [result.MCS.pf;
    result.TDTRA100.pf;
    result.TDTRA150.pf;
    result.TDTRA200.pf;
    result.AMPPT.pf;];
accuracy(:,2) = abs(accuracy(:,1) - result.MCS.pf) ./ result.MCS.pf;

efficiency(:,1) = [result.MCS.funccount;
    result.TDTRA100.funccount;
    result.TDTRA150.funccount;
    result.TDTRA200.funccount;
    result.AMPPT.funccount;];
efficiency(:,2) = [nan;
    100; 150; 200;
    size(result.AMPPT.statistics,1);];
efficiency(:,3) = efficiency(:,1) ./ efficiency(:,2);
table = [accuracy efficiency]; % statistics table


% plot Pf curves
COLOR_MCS = [0, 114, 189] / 255.0;
COLOR_ADMPT = [217, 83, 25] / 255.0;
COLOR_1 = [237, 177, 32] / 255.0;
COLOR_2 = [126, 47, 142] / 255.0;
COLOR_3 = [119, 172, 48] / 255.0;

clf;
plotLine(result.MCS.pf_history, COLOR_MCS, '-');
plotLine(result.TDTRA100.pf_history, COLOR_1, ':v');
plotLine(result.TDTRA150.pf_history, COLOR_2, ':s');
plotLine(result.TDTRA200.pf_history, COLOR_3, ':d');
plotLine(result.AMPPT.pf_history, COLOR_ADMPT,'-o');

scale.figure_width = 8.6;
scale.figure_height = 6.6;
scale.left = 0.14;
scale.bottom = 0.17;
scale.top = 0.08;
scale.right = 0.03;
text.legend_title = {'MCS', 'TDTRA-100', 'TDTRA-150', 'TDTRA-200', 'AMPPT'};
text.x_label = 'Time (month)';
text.y_label = 'Probability of failure';
changeFigureSettings(scale, text);
print('-dtiff','-r300','pf.tif');
end

function changeFigureSettings(scale, text)
set(gca,'fontsize',11, 'fontname','Times New Roman');
set(gcf,'windowstyle','normal');
set(gcf,'unit','centimeters','position',[0 0 scale.figure_width scale.figure_height]);
set(gca,'Position',...
    [scale.left, scale.bottom, ...
    1-scale.left-scale.right, 1-scale.bottom-scale.top]);
h = legend(text.legend_title);
set(h,'FontName','Times New Roman','FontSize',11,'Location','best');
xlabel(text.x_label,'Fontname', 'Times New Roman','FontSize',12);
ylabel(text.y_label,'Fontname', 'Times New Roman','FontSize',12);
box on;
end

function plotLine(data, color, style)
LINE_WIDTH = 1.5;
MARKER_SIZE = 4.0;
hold on;
plot(data(1,:), data(2,:), ...
    style,'LineWidth',LINE_WIDTH, 'MarkerSize', MARKER_SIZE, ...
    'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'Color', color,...
    'MarkerIndices',round(linspace(1,size(data, 2), 10)));
end