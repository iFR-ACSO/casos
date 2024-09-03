clear
close all
clc

addpath(genpath('./pendulum'))

Nmax = 5;
deg  = 2;


%% Run benchmark Nlink pendulum in casos
% disp('Run benchmark in CaSoS')
[bval_array_c, solverTimes_total_c,buildTimes_c] = roaEstNlink_benchCasos(Nmax,deg);

%% Run benchmark Nlink pendulum in spotless
disp('Run benchmark in SPOTless')
[bval_array_sp, solverTimes_total_sp,buildTimes_sp] = roaEstNlink_benchSPOTless(Nmax,deg);


%% plotting
figure('Name','Benchmark N-link Pedulum')
subplot(211)
semilogy((2:Nmax)*2,buildTimes_c,'-o')
hold on
semilogy((2:Nmax)*2,buildTimes_sp,'-+')
% semilogy(2:Nmax,buildTimes_sopt,'-*')
% semilogy(2:Nmax,buildTimes_st,'-^')
xlabel('Number of states')
ylabel('Time [s]')
legend('Ca\Sigmaos','SPOTless')%,'SOSOPT','SOSTOOLS','Location','northwest')
xticks((2:Nmax)*2) 
title('Buildtimes')

subplot(212)
semilogy((2:Nmax)*2,solverTimes_total_c,'-o')
hold on
semilogy((2:Nmax)*2,solverTimes_total_sp,'-+')
% semilogy(2:Nmax,solverTimes_sopt,'-*')
% semilogy(2:Nmax,solverTimes_st,'-^')
xlabel('Number of states')
ylabel('Time [s]')
legend('Ca\Sigmaos','SPOTless')%,'SOSOPT','SOSTOOLS','Location','northwest')
title('Solvertimes')
xticks((2:Nmax)*2)  % Set x-axis ticks to show labels for 2:Nmax