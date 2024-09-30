% ########################################################################
%
% Name:              runBenchmark_III.m
% 
% Short Description: Script that runs the ROA benchmark example for the 
%                    N-link pendulum dynamics described in [1] N-times 
%                    in different toolbox. The mean computation times are 
%                    compared in a bar chart. 
%
% References- [1] Cunis,T., Olucak, J. - CaΣoS: A nonlinear sum-of-squares 
%                 optimization suite, https://arxiv.org/abs/2409.18549
%
% ########################################################################


clear
close all
clc

% add necessary paths
addpath(genpath('./pendulum'))

% add paths to subfolder
addpath('./casos_ex/')
addpath('./sostools_ex/')
addpath('./sosopt_ex/')
addpath('./spotless_ex/')
addpath(genpath('./otherFrameworks'))

% benchmark runs from 2:Nmax
Nmax = 4; 

% degree for polynomial approximation of nonlinear dynamics (currently only
% deg 2 is pre-computed)
deg  = 2;

% we run every N-link test N times, i.e., just set noRuns to specific
% number
noRuns = 5;

%% Run benchmark Nlink pendulum in casos
disp('Run benchmark in CaSoS')
[bval_array_c, solverTimes_total_c,buildTimes_c] = roaEstNlink_benchCasos(Nmax,deg,noRuns);


% we now calculate the mean value of each column i.e. mean over all runs
solverTimes_total_c_mean = mean(solverTimes_total_c,1);
buildTimes_c_mean        = mean(buildTimes_c,1);

%% Run benchmark Nlink pendulum in SOSTOOLS
rmpath(genpath('./otherFrameworks/sosopt'));
addpath(genpath('./otherFrameworks/SOSTOOLS'))
disp('Run benchmark in SOSTOOLS')
[bval_array_st, solverTimes_total_st,buildTimes_st] = roaEstNlink_benchSOSTOOLS(Nmax,deg,noRuns);

% we now calculate the mean value of each column i.e. mean over all runs
solverTimes_total_st_mean = mean(solverTimes_total_st,1);
buildTimes_st_mean        = mean(buildTimes_st,1);


%% Run benchmark Nlink pendulum in spotless
disp('Run benchmark in SPOTless')
[bval_array_sp, solverTimes_total_sp,buildTimes_sp] = roaEstNlink_benchSPOTless(Nmax,deg,noRuns);

% we now calculate the mean value of each column i.e. mean over all runs
solverTimes_total_sp_mean = mean(solverTimes_total_sp,1);
buildTimes_sp_mean        = mean(buildTimes_sp,1);


%% Run benchmark Nlink pendulum in SOSOPT
rmpath(genpath('./otherFrameworks/SOSTOOLS'));
addpath(genpath('./otherFrameworks/sosopt'))
disp('Run benchmark in SOSOPT')
[bval_array_sopt, solverTimes_total_sopt,buildTimes_sopt]= roaEstNlink_benchSosopt(Nmax,deg,noRuns);

% we now calculate the mean value of each column i.e. mean over all runs
solverTimes_total_sopt_mean = mean(solverTimes_total_sopt,1);
buildTimes_sopt_mean        = mean(buildTimes_sopt,1);


save('Benchmark_pendulum_roa.mat')

%% plotting
figure('Name','Benchmark N-link Pedulum')
subplot(211)
semilogy((2:Nmax)*2,buildTimes_c_mean,'-o')
hold on
semilogy((2:Nmax)*2,buildTimes_st_mean,'-^')
semilogy((2:Nmax)*2,buildTimes_sp_mean,'-+')
semilogy((2:Nmax)*2,buildTimes_sopt_mean,'-*')
xlabel('Number of states')
ylabel('Parsing time [s]')
legend('Ca\Sigmaos','SOSTOOLS','SPOTless','SOSOPT','Location','northwest')
xticks((2:Nmax)*2) 


subplot(212)
semilogy((2:Nmax)*2,solverTimes_total_c_mean,'-o')
hold on
semilogy((2:Nmax)*2,solverTimes_total_st_mean,'-^')
semilogy((2:Nmax)*2,solverTimes_total_sp_mean,'-+')
semilogy((2:Nmax)*2,solverTimes_total_sopt_mean,'-*')

xlabel('Number of states')
ylabel('Solver time [s]')
legend('Ca\Sigmaos','SOSTOOLS','SPOTless','SOSOPT','Location','northwest')
xticks((2:Nmax)*2) 

cleanfigure()
matlab2tikz('benchmark_pendulum_roa.tex','width','\figW','height','\figH');