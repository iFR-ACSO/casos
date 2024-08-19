clc
clear
close all


% maximum number of states
Nmax = 2;

% degree of Taylor-Approximation; currently only pre-computed for deg = 2!
deg  = 2;

% solve N-link Pendulum ROA with CaSoS
[solverTimes_c,buildTimes_c,gval_arr_c] = roaEstNlinkPend_benchVsIterCasos(Nmax,deg);


% solve N-link Pendulum ROA with SPOTless
[solverTimes_sopt,buildTimes_sopt,gval_arr_sopt]= roaEstNlinkPend_benchVsIterSOSOPT(Nmax,deg);


[solverTimes_sp,buildTimes_sp,gval_arr_sp]= roaEstNlinkPend_benchVsIterSPOTless(Nmax,deg);
% ------------------------------------------------------------------------
% some functions in SOSTOOLS and SOSOPT have same names; since both use
% multipoly it might not work; remove one off them an make sure to add the
% other one
% ------------------------------------------------------------------------
rmpath(genpath('C:\Users\ac133867\Documents\MATLAB\sosopt'));
addpath(genpath('C:\Users\ac133867\Documents\MATLAB\SOSTOOLS'))
[solverTimes_st,buildTimes_st,gval_array_st]= roaEstNlinkPend_benchSOSTOOLS(Nmax,deg);

% solve N-link Pendulum ROA with SOSOPT
% ------------------------------------------------------------------------
% some functions in SOSTOOLS and SOSOPT have same names; since both use
% multipoly it might not work; remove one off them an make sure to add the
% other one
% ------------------------------------------------------------------------
rmpath(genpath('C:\Users\ac133867\Documents\MATLAB\SOSTOOLS'));
addpath(genpath('C:\Users\ac133867\Documents\MATLAB\sosopt'))
% [solverTimes_sopt,buildTimes_sopt,gval_arr_sopt]= roaEstNlinkPend_benchSOSOPT(Nmax,deg);

%% plotting
figure('Name','Benchmark N-link Pedulum')
subplot(211)
semilogy(2:Nmax,buildTimes_c,'-o')
hold on
semilogy(2:Nmax,buildTimes_sp,'-+')
semilogy(2:Nmax,buildTimes_sopt,'-*')
semilogy(2:Nmax,buildTimes_st,'-^')
xlabel('Number of states')
ylabel('Time [s]')
legend('Ca\Sigmaos','SPOTless','SOSOPT','SOSTOOLS','Location','northwest')
title('Buildtimes')
subplot(212)
semilogy(2:Nmax,solverTimes_c,'-o')
hold on
semilogy(2:Nmax,solverTimes_sp,'-+')
semilogy(2:Nmax,solverTimes_sopt,'-*')
semilogy(2:Nmax,solverTimes_st,'-^')
xlabel('Number of states')
ylabel('Time [s]')
legend('Ca\Sigmaos','SPOTless','SOSOPT','SOSTOOLS','Location','northwest')
title('Solvertimes')