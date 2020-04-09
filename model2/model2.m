function model2()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% problem definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
problem.performanceFunc = @CantileverBeam;
problem.Ts = 0;
problem.Te = 5;
problem.variable_table = {
    % distribution      mean        std             autocorrelationCoef
    
    % d, h, R0, F2, P
    'normal',           0.042,      0.042*0.0119,   [];
    'normal',           0.005,      0.005*0.02,     [];
    'normal',           560e6,      560e5,          [];
    'normal',           1.8e3,      1.8e2,          [];
    'gumbel',           1.0e3,      1.0e2,          [];
    
    % F1 T(t)
    'stationary',       1.8e3,      1.8e2,          @(t1,t2) exp(-abs(t1-t2)/4);
    'stationary',       1.9e3,      1.9e2,          @(t1,t2) exp(-4*(t1-t2).^2);
    };

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
option.discretization_node_num = 100;
result.AMPPT = AMPPT( problem, option );

% TDTRA
option.time_node_num = 10; result.TDTRA10 = TDTRA( problem, option );
option.time_node_num = 20; result.TDTRA20 = TDTRA( problem, option );
option.time_node_num = 30; result.TDTRA30 = TDTRA( problem, option );
save('result.mat', 'result');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% analyze the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analyzeResults()
result = [];
load('result.mat');

accuracy = [result.MCS.pf;
    result.TDTRA10.pf;
    result.TDTRA20.pf;
    result.TDTRA30.pf;
    result.AMPPT.pf;];
accuracy(:,2) = abs(accuracy(:,1) - result.MCS.pf) ./ result.MCS.pf;

efficiency(:,1) = [result.MCS.funccount;
    result.TDTRA10.funccount;
    result.TDTRA20.funccount;
    result.TDTRA30.funccount;
    result.AMPPT.funccount;];
efficiency(:,2) = [nan;
    10; 20; 30;
    5;];
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
plotLine(result.TDTRA10.pf_history, COLOR_1, ':v');
plotLine(result.TDTRA20.pf_history, COLOR_2, ':s');
plotLine(result.TDTRA30.pf_history, COLOR_3, ':d');
plotLine(result.AMPPT.pf_history, COLOR_ADMPT,'-o');

scale.figure_width = 8.6;
scale.figure_height = 6.6;
scale.left = 0.18;
scale.bottom = 0.17;
scale.margin = 0.02;
text.legend_title = {'MCS', 'TDTRA-10', 'TDTRA-20', 'TDTRA-30', 'AMPPT'};
text.x_label = 'Time (year)';
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
    1-scale.left-scale.margin, 1-scale.bottom-scale.margin]);

h = legend(text.legend_title);
set(h,'FontName','Times New Roman','FontSize',12,'Location','best');

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
