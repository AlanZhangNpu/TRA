function model1()
rng('shuffle');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% problem definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
problem.performanceFunc = @SteelBeam;
problem.Ts = 0;
problem.Te = 30;
problem.variable_table = {
    % distribution      mean        std         autocorrelationCoef
    % w0, h0, sigma
    'lognormal',        0.20,       0.01,       [];
    'lognormal',        0.04,       4e-3,       [];
    'lognormal',        2.4e8,      2.4e7,      [];
    % F(t)
    'stationary',       3500,       700,        @(t1, t2) exp(-(t1-t2).^2);
    };

solve(problem);
analyzeResults();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% solve the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function solve(problem)
% MCS
MCS_par.n_MCS = 2e6;
MCS_par.n_time_instant = 500;
result.MCS = MCS(problem, MCS_par);

% ADMPT
option.discretization_node_num = 70;
result.AMPPT = AMPPT( problem, option );

% TDRA
option.time_node_num = 30;  result.TDRA30 = TDTRA( problem, option );
option.time_node_num = 40;  result.TDRA40 = TDTRA( problem, option );
option.time_node_num = 50;  result.TDRA50 = TDTRA( problem, option );
save('result.mat', 'result');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% analyze the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function analyzeResults()
result = [];
load('result.mat');

accuracy = [result.MCS.pf;
    result.TDTRA30.pf;
    result.TDTRA40.pf;
    result.TDTRA50.pf;
    result.AMPPT.pf;];
accuracy(:,2) = abs(accuracy(:,1) - result.MCS.pf) ./ result.MCS.pf;

efficiency(:,1) = [result.MCS.funccount;
    result.TDTRA30.funccount;
    result.TDTRA40.funccount;
    result.TDTRA50.funccount;
    result.AMPPT.funccount;];
efficiency(:,2) = [nan; 30; 40; 50; 5;];
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
plotLine(result.TDTRA30.pf_history, COLOR_1, ':v');
plotLine(result.TDTRA40.pf_history, COLOR_2, ':s');
plotLine(result.TDTRA50.pf_history, COLOR_3, ':d');
plotLine(result.AMPPT.pf_history, COLOR_ADMPT,'-o');

scale.figure_width = 8.6;
scale.figure_height = 6.6;
scale.left = 0.14;
scale.bottom = 0.17;
scale.top = 0.08;
scale.right = 0.03;
text.legend_title = {'MCS', 'TDTRA-30', 'TDTRA-40', 'TDTRA-50', 'AMPPT'};
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
    1-scale.left-scale.right, 1-scale.bottom-scale.top]);
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