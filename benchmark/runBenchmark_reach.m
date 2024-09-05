% script that runs the GTM 4D ROA problem with different toolboxes

clc
clear
close all

% add paths to subfolder
addpath('./casos_ex/')
addpath('./yalmip_ex/')
addpath('./sostools_ex/')
addpath('./sosopt_ex/')
addpath('./spotless_ex/')

%% run GTM 4D example in casos
disp('Run benchmark test for CaSoS')
[gval_c_GTM,solverTime_total_c_GTM,buildTime_c_GTM]= reachEstGTM_benchCasos();

% %% run GTM 4D example in sosopt
% disp('Run benchmark test for SOSOPT')
% rmpath(genpath('C:\Users\ac133867\Documents\MATLAB\SOSTOOLS'));
% addpath(genpath('C:\Users\ac133867\Documents\MATLAB\sosopt'))
[gval_sopt_GTM,solverTime_total_sopt_GTM,buildTime_sopt_GTM] = reachEstGTM_benchSosopt();

%% run GTM 4D example in spotless
disp('Run benchmark test for SPOTless')
[gval_sp_GTM,solverTime_sp_GTM,buildTime_sp_GTM]= reachEstGTM_benchSPOTless();

%% run GTM 4D example in SOSTOOLS
disp('Run benchmark test for SOSTOOLS')
% rmpath(genpath('C:\Users\ac133867\Documents\MATLAB\sosopt'));
% addpath(genpath('C:\Users\ac133867\Documents\MATLAB\SOSTOOLS'))
[gval_st_GTM,solverTime_st_GTM,buildTime_st_GTM]= reachEstGTM_benchSOSTOOLS();

%% plot result
solvers = {'Ca\SigmaoS' ,'SOSOPT','SPOTless','SOSTOOLS'};

buildTimes  = [buildTime_c_GTM,buildTime_sopt_GTM,buildTime_sp_GTM,buildTime_st_GTM];                                          % Time spent building for each solver
solverTimes = [solverTime_total_c_GTM,solverTime_total_sopt_GTM, solverTime_sp_GTM,solverTime_st_GTM]; % Time spend in low-level solver 

% Combine the build and solve times into a matrix
timeData = [buildTimes; solverTimes]';

% Create the stacked bar chart
figure('Name','Comparison solve and build time');
bar(timeData, ... % actual data
    0.4, ...      % width of bars
    'stacked');   % type

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

axis([0.5 length(solvers)+0.5 0 max(totalTimes)+0.1*max(totalTimes)])
% Add a legend to indicate what each part of the stack represents
legend('Parsing/Build Time', 'Solver Time','Location', 'northwest');