% ########################################################################
%
% Name:              runBenchmark_pendulum.m
% 
% Short Description: Script that runs the ROA benchmark example for the 
%                    N-link pendulum dynamics described in [1] N-times 
%                    in different toolbox. The mean computation times are 
%                    compared in a bar chart. 
%
% References-
%
% ########################################################################


clear
close all
clc

% add necessary paths
addpath(genpath('./pendulum'))

% add paths to subfolder
addpath('./casos_ex/')
addpath('./yalmip_ex/')
addpath('./sostools_ex/')
addpath('./sosopt_ex/')
addpath('./spotless_ex/')
addpath(genpath('./otherFrameworks'))

% benchmark runs from 2:Nmax
Nmax = 4;

% degree for polynomial approximation of nonlinear dynamics
deg  = 2;

%% Run benchmark Nlink pendulum in casos
disp('Run benchmark in CaSoS')
[bval_array_c, solverTimes_total_c,buildTimes_c] = roaEstNlink_benchCasos(Nmax,deg);

%% Run benchmark Nlink pendulum in spotless
disp('Run benchmark in SPOTless')
[bval_array_sp, solverTimes_total_sp,buildTimes_sp] = roaEstNlink_benchSPOTless(Nmax,deg);

%% Run benchmark Nlink pendulum in SOSTOOLS
rmpath(genpath('./otherFrameworks/sosopt'));
addpath(genpath('./otherFrameworks/SOSTOOLS'))
disp('Run benchmark in SOSTOOLS')
[bval_array_st, solverTimes_total_st,buildTimes_st] = roaEstNlink_benchSOSTOOLS(Nmax,deg);

%% Run benchmark Nlink pendulum in SOSOPT
rmpath(genpath('./otherFrameworks/SOSTOOLS'));
addpath(genpath('./otherFrameworks/sosopt'))
disp('Run benchmark in SOSOPT')
[bval_array_sopt, solverTimes_total_sopt,buildTimes_sopt]= roaEstNlink_benchSosopt(Nmax,deg);


save('Benchmark_pendulum_roa.mat')
%% plotting
figure('Name','Benchmark N-link Pedulum')
subplot(211)
semilogy((2:Nmax)*2,buildTimes_c,'-o')
hold on
semilogy((2:Nmax)*2,buildTimes_sp,'-+')
semilogy((2:Nmax)*2,buildTimes_sopt,'-*')
semilogy((2:Nmax)*2,buildTimes_st,'-^')
xlabel('Number of states')
ylabel('Time [s]')
legend('Ca\Sigmaos','SPOTless','SOSOPT','SOSTOOLS','Location','northwest')
xticks((2:Nmax)*2) 
% title('Buildtimes')

subplot(212)
semilogy((2:Nmax)*2,solverTimes_total_c,'-o')
hold on
semilogy((2:Nmax)*2,solverTimes_total_sp,'-+')
semilogy((2:Nmax)*2,solverTimes_total_sopt,'-*')
semilogy((2:Nmax)*2,solverTimes_total_st,'-^')
xlabel('Number of states')
ylabel('Time [s]')
legend('Ca\Sigmaos','SPOTless','SOSOPT','SOSTOOLS','Location','northwest')%,'SOSOPT','SOSTOOLS','Location','northwest')
% title('Solvertimes')
xticks((2:Nmax)*2)  % Set x-axis ticks to show labels for 2:Nmax

cleanfigure()
matlab2tikz('benchmark_pendulum_roa.tex','width','\figW','height','\figH');