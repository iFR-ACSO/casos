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

%% Short description benchmark
% make five runs to get a better mean value of the computation times; 
% keep the last gamma and beta value 

% pre-allocate
solverTime_total_c_GTM_arr = zeros(5,1);
buildTime_c_GTM_arr        = zeros(5,1);
callTime_c_GTM_arr        = zeros(5,1);

solverTime_total_st_GTM_arr = zeros(5,1);
buildTime_st_GTM_arr        = zeros(5,1);

solverTime_total_st2_GTM_arr = zeros(5,1);
buildTime_st2_GTM_arr        = zeros(5,1);

solverTime_total_sp_GTM_arr = zeros(5,1);
buildTime_sp_GTM_arr        = zeros(5,1);

solverTime_total_sopt_GTM_arr = zeros(5,1);
buildTime_sopt_GTM_arr        = zeros(5,1);

solverTime_total_y_GTM_arr = zeros(5,1);
buildTime_y_GTM_arr        = zeros(5,1);


for j = 1:1

%% run GTM 4D example in casos
disp('Run benchmark test for CaSoS')
[gval_c_GTM,bval_c_GTM,solverTime_total_c_GTM,buildTime_c_GTM,callTime_c_GTM] = roaEstGTM_benchCasos();

solverTime_total_c_GTM_arr(j)  = solverTime_total_c_GTM;
buildTime_c_GTM_arr(j)         = buildTime_c_GTM;
callTime_c_GTM_arr(j)          = callTime_c_GTM;

% %% run GTM 4D example in SOSTOOLS with dpvar and pvar
% disp('Run benchmark test for SOSTOOLS with dpvar')
% % ------------------------------------------------------------------------
% % some functions in SOSTOOLS and SOSOPT have same names; since both use
% % multipoly it might not work; remove one off them an make sure to add the
% % other one
% % ------------------------------------------------------------------------
% rmpath(genpath('./otherFrameworks/sosopt'));
% addpath(genpath('./otherFrameworks/SOSTOOLS'))
% [gval_st_GTM,bval_st_GTM,solverTime_total_st_GTM,buildTime_st_GTM]= roaEstGTM_benchSOSTOOLS();
% 
% solverTime_total_st_GTM_arr(j)  = solverTime_total_st_GTM;
% buildTime_st_GTM_arr(j)         = buildTime_st_GTM;
% 
% % run GTM 4D example in SOSTOOLS with pvar
% disp('Run benchmark test for SOSTOOLS with pvar')
% [gval_st2_GTM,bval_st2_GTM,solverTime_total_st2_GTM,buildTime_st2_GTM]= roaEstGTM_benchSOSTOOLS2();
% 
% solverTime_total_st2_GTM_arr(j)  = solverTime_total_st2_GTM;
% buildTime_st2_GTM_arr(j)         = buildTime_st2_GTM;
% 
% %% run GTM 4D example in sosopt using gsosopt
% disp('Run benchmark test for SOSOPT using GSOSOPT')
% ------------------------------------------------------------------------
% some functions in SOSTOOLS and SOSOPT have same names; since both use
% multipoly it might not work; remove one off them an make sure to add the
% other one
% % ------------------------------------------------------------------------
% rmpath(genpath('./otherFrameworks/SOSTOOLS'));
% addpath(genpath('./otherFrameworks/sosopt'))
% [gval_sopt_GTM,bval_sopt_GTM,solverTime_total_sopt_GTM,buildTime_sopt_GTM]= roaEstGTM_benchSosopt();
% 
% solverTime_total_sopt_GTM_arr(j)  = solverTime_total_sopt_GTM;
% buildTime_sopt_GTM_arr(j)         = buildTime_sopt_GTM;
% 
% %% run GTM 4D example in spotless
% disp('Run benchmark test for SPOTless')
% % run GTM 4D example in yalmip with default options
% [gval_sp_GTM,bval_sp_GTM,solverTime_total_sp_GTM,buildTime_sp_GTM] = roaEstGTM_benchSPOTless();
% 
% solverTime_total_sp_GTM_arr(j)  = solverTime_total_sp_GTM;
% buildTime_sp_GTM_arr(j)         = buildTime_sp_GTM;
% 
% %% run GTM 4D example in YALMIP
% disp('Run benchmark test for Yalmip')
% % run GTM 4D example in yalmip with default options
% [gval_y_GTM,bval_y_GTM,solverTime_total_y_GTM,buildTime_y_GTM] = roaEstGTM_benchYALMIP();
% 
% solverTime_total_y_GTM_arr(j) = solverTime_total_y_GTM;
% buildTime_y_GTM_arr(j)        = buildTime_y_GTM;

end

% from the five runs, take the mean value
solverTime_total_c_GTM     = mean(solverTime_total_c_GTM_arr);
solverTime_total_st_GTM    = mean(solverTime_total_st_GTM_arr);
solverTime_total_sp_GTM    = mean(solverTime_total_sp_GTM_arr);
solverTime_total_sopt_GTM  = mean(solverTime_total_sopt_GTM_arr);
solverTime_total_y_GTM     = mean(solverTime_total_y_GTM_arr);

buildTime_c_GTM     = mean(buildTime_c_GTM_arr);
buildTime_st_GTM    = mean(buildTime_st_GTM_arr);
buildTime_sp_GTM    = mean(buildTime_sp_GTM_arr);
buildTime_sopt_GTM  = mean(buildTime_sopt_GTM_arr);
buildTime_y_GTM     = mean(buildTime_y_GTM_arr);

% save latest benchmark
save('Benchmark_gtm_roa.mat')

%% plot result
solvers = {'Ca\SigmaoS', 'SOSTOOLS (dpvar)', 'SOSTOOLS (pvar)' ,'SOSOPT','SPOTless' ,'YALMIP'};

buildTimes  = [buildTime_c_GTM,buildTime_st_GTM,buildTime_st2_GTM, buildTime_sopt_GTM,buildTime_sp_GTM, buildTime_y_GTM];                                          % Time spent building for each solver
solverTimes = [solverTime_total_c_GTM,solverTime_total_st_GTM,solverTime_total_st2_GTM, solverTime_total_sopt_GTM,solverTime_total_sp_GTM,solverTime_total_y_GTM]; % Time spend in low-level solver 

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

axis([0.5 length(solvers)+0.5 0 1600])
% Add a legend to indicate what each part of the stack represents
legend('Parsing/Build Time', 'Solver Time','Location', 'northwest');

cleanfigure()
matlab2tikz('benchmark_gtm_roa.tex','width','\figW','height','\figH');