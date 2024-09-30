% ########################################################################
%
% Name:              runBenchmark_reach.m
% 
% Short Description: Script that runs the reachability benchmark example 
%                    described in [1] N-times indifferent toolbox. The mean 
%                    computation times are compared in a bar chart. The
%                    implemented SOS problem and derivation can be found in
%                    [2]. The original implemention of [2] can be found
%                    here:https://github.com/heyinUCB/Backward-Reachability-
%                    Analysis-and-Control-Synthesis. See comments in the
%                    SOSOPT benchmark implementation.
%
% References- [1] Cunis,T., Olucak, J. - CaΣoS: A nonlinear sum-of-squares 
%                 optimization suite
%             [2] Yin, H. et al- Backward Reachability for Polynomial 
%                 Systems on a Finite Horizon, IEEE Transactions on 
%                 Automatic Control, 2021, doi: 10.1109/TAC.2021.3056611
%
% ########################################################################

clc
clear
close all

% add paths to subfolder
addpath('./casos_ex/')
addpath('./yalmip_ex/')
addpath('./sostools_ex/')
addpath('./sosopt_ex/')
addpath('./spotless_ex/')
addpath(genpath('./otherFrameworks'))


%% Short description benchmark
% make  n number of runs (noRuns) to get a better mean value of the computation times;
% keep the last gamma value and the last storage function for later volume
% computation

noRuns = 5;

% pre-allocate
solverTime_total_c_GTM_arr = zeros(noRuns,1);
buildTime_c_GTM_arr        = zeros(noRuns,1);

solverTime_total_st_GTM_arr = zeros(noRuns,1);
buildTime_st_GTM_arr        = zeros(noRuns,1);

solverTime_total_sp_GTM_arr = zeros(noRuns,1);
buildTime_sp_GTM_arr        = zeros(noRuns,1);

solverTime_total_sopt_GTM_arr = zeros(noRuns,1);
buildTime_sopt_GTM_arr        = zeros(noRuns,1);

for j = 1:noRuns

disp(['Round number: ' num2str(j)])
%% run GTM 4D example in casos
disp('Run benchmark test for CaSoS')
[gval_c_GTM,solverTime_total_c_GTM,buildTime_c_GTM,V_c]= reachEstGTM_benchCasos();

solverTime_total_c_GTM_arr(j)  = solverTime_total_c_GTM;
buildTime_c_GTM_arr(j)         = buildTime_c_GTM;

%% run GTM 4D example in SOSTOOLS
disp('Run benchmark test for SOSTOOLS')
% ------------------------------------------------------------------------
% some functions in SOSTOOLS and SOSOPT have same names; since both use
% multipoly it might not work; remove one off them an make sure to add the
% other one
% ------------------------------------------------------------------------
rmpath(genpath('./otherFrameworks/sosopt'));
addpath(genpath('./otherFrameworks/SOSTOOLS'))
[gval_st_GTM,solverTime_st_GTM,buildTime_st_GTM,V_st]= reachEstGTM_benchSOSTOOLS();

solverTime_total_st_GTM_arr(j)  = solverTime_st_GTM;
buildTime_st_GTM_arr(j)         = buildTime_st_GTM;

%% run GTM 4D example in spotless
disp('Run benchmark test for SPOTless')
[gval_sp_GTM,solverTime_sp_GTM,buildTime_sp_GTM,V_sp]= reachEstGTM_benchSPOTless();

solverTime_total_sp_GTM_arr(j)  = solverTime_sp_GTM;
buildTime_sp_GTM_arr(j)         = buildTime_sp_GTM;

%% run GTM 4D example in sosopt
disp('Run benchmark test for SOSOPT')
% ------------------------------------------------------------------------
% some functions in SOSTOOLS and SOSOPT have same names; since both use
% multipoly it might not work; remove one off them an make sure to add the
% other one
% ------------------------------------------------------------------------
rmpath(genpath('./otherFrameworks/SOSTOOLS'));
addpath(genpath('./otherFrameworks/sosopt'))
[gval_sopt_GTM,solverTime_sopt_GTM,buildTime_sopt_GTM,V_sopt] = reachEstGTM_benchSosopt();

solverTime_total_sopt_GTM_arr(j)  = solverTime_sopt_GTM;
buildTime_sopt_GTM_arr(j)         = buildTime_sopt_GTM;

end

% save latest benchmark
save('Benchmark_gtm_reach.mat')

%% get mean computation times
% from the five runs, take the mean value
buildTime_c_GTM            = mean(buildTime_c_GTM_arr);
solverTime_total_c_GTM     = mean(solverTime_total_c_GTM_arr);
total_c                    = buildTime_c_GTM + solverTime_total_c_GTM;

buildTime_st_GTM           = mean(buildTime_st_GTM_arr);
solverTime_total_st_GTM    = mean(solverTime_total_st_GTM_arr);
total_st                   = buildTime_st_GTM + solverTime_total_st_GTM;

buildTime_sp_GTM           = mean(buildTime_sp_GTM_arr);
solverTime_total_sp_GTM    = mean(solverTime_total_sp_GTM_arr);
total_sp                   = buildTime_sp_GTM + solverTime_total_sp_GTM;

buildTime_sopt_GTM         = mean(buildTime_sopt_GTM_arr);
solverTime_total_sopt_GTM  = mean(solverTime_total_sopt_GTM_arr);
total_sopt                 = buildTime_sopt_GTM + solverTime_total_sopt_GTM;


%% plot result
% solvers = {'Ca\SigmaoS' ,'SOSTOOLS' , 'SPOTless','SOSOPT'};
solvers = {'A' ,'B' , 'C','D'};
% Time spent building/parsing for each solver
buildTimes  = [buildTime_c_GTM, buildTime_st_GTM, buildTime_sp_GTM, buildTime_sopt_GTM];           
% Time spend in low-level solver 
solverTimes = [solverTime_total_c_GTM, solverTime_total_st_GTM, solverTime_total_sp_GTM, solverTime_total_sopt_GTM]; 

% Combine the build and solve times into a matrix
timeData = [buildTimes; solverTimes]';

% Create the stacked bar chart
figure('Name','Comparison solve and build time');
bar(timeData, ... % actual data
    0.4, ...      % width of bars
    'stacked');   % type

% Customize the x-axis with solver names
set(gca, 'XTickLabel', solvers);

% Add labels 
ylabel('Time (seconds)');

% Display the figure
grid on;

% Calculate the total times for each solver (buildTime + solverTime)
totalTimes = sum(timeData, 2);

axis([0.5 length(solvers)+0.5 0 max(totalTimes)+0.1*max(totalTimes)])
% Add a legend to indicate what each part of the stack represents
legend('Parsing/Build Time', 'Solver Time','Location', 'northwest');

%% make tikz figure (matlab2tikz must be on the matlab path)
% cleanfigure()
% matlab2tikz('benchmark_gtm_reach.tex','width','\figW','height','\figH');

%% plotting of sublevel sets; here slice of scaled x1-x2
import casos.toolboxes.to_multipoly
pvar x_1 x_2 x_3 x_4 t % for plotting
pvar x1 x2 x3 x4       % for plotting

% casos
Vc_indet = V_c.indeterminates;
V0_c = subs(V_c,Vc_indet(1),0); % V(0,x)

% we use multipoly function to compute the volume also for CaSoS
V_c_mult =  to_multipoly(V0_c - gval_c_GTM);

figure(3)
pcontour(subs(V_c_mult,[x_3;x_4],[0;0]),0,[-1 1 -1 1],'r')
hold on
% SOSTOOLS and SOSOPT both use already multipoly
pcontour(subs(subs(V_sopt,[x_3;x_4],[0;0]),t,0)-gval_sopt_GTM,0,[-1 1 -1 1],'b')
pcontour(subs(subs(V_st,[x3;x4],[0;0]),t,0)-gval_st_GTM,0,[-1 1 -1 1],'g')
legend('Ca')

% spotless needs a custom solution
x = msspoly ('x' , 4 ) ; % for plotting
t = msspoly ('t' , 1 ) ;
V0_sp = subs(V_sp,t,0);

% plot spotless solution
domain = [-5 5 -5 5]./5;
xg      = linspace(domain(1),domain(2),100);
yg      = linspace(domain(3),domain(4),100);
[xg,yg] = meshgrid(xg,yg);

V0_sp = subs(V0_sp,x(3:4),[0;0]);
pgrid = double(msubs(V0_sp,x(1:2),[xg(:)'; yg(:)']));
pgrid = reshape(pgrid,size(xg));
contour(xg,yg,double(pgrid),[0,0],'k');
grid minor
legend('Ca\SigmaoS' ,'SOSTOOLS' ,'SOSOPT' , 'SPOTless')

%% estimate volume via Monte-Carlo-Sampling (!!! result might differ due to randomness!!!)
% compute volume; all four use a unit box as bounding region

[vol_sp,std_sp] = pvolume_spot(V0_sp,x,gval_sp_GTM);  % spot

% sostools and sosopt use both multipoly; we transformed CaSoS to multipoly
[vol_c,std_c]       = pvolume(V_c_mult,0); % casos
pvar t
[vol_st,std_st]     = pvolume(subs(V_st,t,0)-gval_st_GTM,0); % sostools
[vol_sopt,std_sopt] = pvolume(subs(V_sopt,t,0)-gval_sopt_GTM,0); % sosopt
