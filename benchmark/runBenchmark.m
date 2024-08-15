
clc
clear
close all

addpath('./casos_ex/')
addpath('./yalmip_ex/')
addpath('./sostools_ex/')
addpath('./sosopt_ex/')
addpath('./spotless_ex/')
% run GTM 4D example in casos
disp('Run benchmark test for CaSoS')
[gval_c_GTM,bval_c_GTM,solveTime_total_c_GTM,solverTime_total_c_GTM,buildTime_c_GTM] = roaEstGTM_benchCasos();

disp('Run benchmark test for SOSTOOLS with dpvar')
rmpath(genpath('C:\Users\ac133867\Documents\MATLAB\sosopt'));
addpath(genpath('C:\Users\ac133867\Documents\MATLAB\SOSTOOLS'))
[gval_st_GTM,bval_st_GTM,solveTime_total_st_GTM,solverTime_total_st_GTM,buildTime_st_GTM]= roaEstGTM_benchSOSTOOLS();

disp('Run benchmark test for SOSTOOLS with pvar')
[gval_st2_GTM,bval_st2_GTM,solveTime_total_st2_GTM,solverTime_total_st2_GTM,buildTime_st2_GTM]= roaEstGTM_benchSOSTOOLS2();

disp('Run benchmark test for SOSOPT using GSOSOPT')
rmpath(genpath('C:\Users\ac133867\Documents\MATLAB\SOSTOOLS'));
addpath(genpath('C:\Users\ac133867\Documents\MATLAB\sosopt'))
[gval_sopt_GTM,bval_sopt_GTM,solveTime_total_sopt_GTM,solverTime_total_sopt_GTM,buildTime_sopt_GTM]= roaEstGTM_benchSosopt();


disp('Run benchmark test for SPOTless')
% run GTM 4D example in yalmip with default options
[gval_sp_GTM,bval_sp_GTM,solveTime_total_sp_GTM,solverTime_total_sp_GTM,buildTime_sp_GTM] = roaEstGTM_benchSPOTless();

disp('Run benchmark test for Yalmip')
% run GTM 4D example in yalmip with default options
defaultOpts = 0;
[gval_y_GTM,bval_y_GTM,solveTime_total_y_GTM,solverTime_total_y_GTM,buildTime_y_GTM] = roaEstGTM_benchYALMIP(defaultOpts);



%% plot result
solvers = {'Ca\SigmaoS', 'SOSTOOLS (dpvar)', 'SOSTOOLS (pvar)' ,'SOSOPT','SPOTless' ,'YALMIP'};

buildTimes  = [buildTime_c_GTM,buildTime_st_GTM,buildTime_st2_GTM, buildTime_sopt_GTM,buildTime_sp_GTM, buildTime_y_GTM];               % Time spent building for each solver
solveTimes  = [solveTime_total_c_GTM, solveTime_total_st_GTM,solveTime_total_st2_GTM, solveTime_total_sopt_GTM,solveTime_total_sp_GTM, solveTime_total_y_GTM];   % Time spent solving for each quasi-convex problem
solverTimes = [solverTime_total_c_GTM,solverTime_total_st_GTM,solverTime_total_st2_GTM, solverTime_total_sopt_GTM,solverTime_total_sp_GTM,solverTime_total_y_GTM]; % Time spend in low-level solver 

% Combine the build and solve times into a matrix
timeData = [buildTimes; solverTimes]';

% Create the stacked bar chart
figure('Name','Comparison solve and build time');
bar(timeData,0.4, 'stacked');

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

axis([0.5 length(solvers)+0.5 0 1300])
% Add a legend to indicate what each part of the stack represents
legend('Parsing/Build Time', 'Solver Time','Location', 'northwest');

cleanfigure()
matlab2tikz('benchmark_gtm_roa.tex','width','\figW','height','\figH');