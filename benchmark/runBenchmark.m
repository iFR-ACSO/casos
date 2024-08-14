
clc
clear
close all

addpath('./casos_ex/')
addpath('./yalmip_ex/')

% run GTM 4D example in casos
[gval_c_GTM,bval_c_GTM,solveTime_total_c_GTM,solverTime_total_c_GTM,buildTime_c_GTM] = roaEstGTM_benchCasos();

% run GTM 4D example in yalmip with default options
defaultOpts = 0;
[gval_y_GTM,bval_y_GTM,solveTime_total_y_GTM,solverTime_total_y_GTM,buildTime_y_GTM] = roaEstGTM_benchYALMIP(defaultOpts);


%% plot result
solvers = {'Ca\SigmaoS', 'YALMIP'};

buildTimes  = [buildTime_c_GTM, buildTime_y_GTM];               % Time spent building for each solver
solveTimes  = [solveTime_total_c_GTM, solveTime_total_y_GTM];   % Time spent solving for each quasi-convex problem
solverTimes = [solverTime_total_c_GTM, solverTime_total_y_GTM]; % Time spend in low-level solver 

% Combine the build and solve times into a matrix
timeData = [buildTimes; solverTimes]';

% Create the stacked bar chart
figure('Name','Comparison solve and build time');
bar(timeData, 'stacked');

% Customize the x-axis with solver names
set(gca, 'XTickLabel', solvers);

% Add labels and title
ylabel('Time (seconds)');

% Display the figure
grid on;

% Calculate the total times for each solver (buildTime + solverTime)
totalTimes = sum(timeData, 2);

% Add the total times as text labels at the top of each stack
x = 1:length(solvers); % X positions for each bar

% Loop through each bar to place the text
for i = 1:length(x)
    % The position is above the top of each stack (build + solve)
    text(x(i), totalTimes(i), sprintf('%.2f s', totalTimes(i)), ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end


% Add a legend to indicate what each part of the stack represents
legend('Parsing/Build Time', 'Solver Time','Location', 'northwest');


