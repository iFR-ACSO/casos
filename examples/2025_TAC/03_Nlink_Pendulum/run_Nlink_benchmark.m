clc
clear


deg    = 2;
noRuns = 1;
Nmax   = 10; 

startTotal = tic;

% run sequential SOS
[fval_array_Seq, solverTimes_total_Seq,buildTimes_Seq,iteration_array_Seq,solStatus_array_Seq,solverStats_Seq] = bench_Nlink_Seq(noRuns,deg,Nmax);

% run coordinate descent
[fval_array_CD, solverTimes_total_CD,buildTimes_CD,iteration_array_CD,solverStats_CD] = bench_Nlink_CD(noRuns,deg,Nmax);


% run sequential SOS bad initial
[fval_array_Seq_badInit, solverTimes_total_Seq_badInit,buildTimes_Seq_badInit,iteration_array_Seq_badInit,solStatus_array_Seq_badInit,solverStats_Seq_badInit] = bench_Nlink_Seq_badInit(noRuns,deg,Nmax);

totalBenchTime = toc(startTotal);


save('Nlink_alsoBadinint_2802.mat')

%% evaluation
runs = 2:Nmax;

% plot solve times
figure(1)
subplot(211)
semilogy(runs*2,solverTimes_total_Seq,'-o')
hold on
semilogy(runs*2,solverTimes_total_Seq_badInit,'-or')
semilogy(runs*2,solverTimes_total_CD,'-^')

% xlabel('Number of states')
ylabel('Solver time [s]')
legend('sequential','sequential (bad initial guess)','coordinate-descent','Location','northwest')
xticks((2:Nmax)*2) 

% plot iterations (no subiterations of feasibiltity restoration)
subplot(212)
plot(runs*2,iteration_array_Seq,'-o')
hold on
plot(runs*2,iteration_array_Seq_badInit,'-or')
plot(runs*2,iteration_array_CD,'-^')

xlabel('Number of states')
ylabel('Iterations [-]')
legend('sequential','sequential (bad initial guess)','coordinate-descent','Location','northwest')
xticks((2:Nmax)*2) 


cleanfigure();
matlab2tikz('NlinkPendulum.tex','width','\figW','height','\figH')

% time per iteration
figure(3)
plot(runs*2,solverTimes_total_Seq./iteration_array_Seq,'-o')
hold on
plot(runs*2,solverTimes_total_Seq_badInit./iteration_array_Seq_badInit,'-o')
plot(runs*2,solverTimes_total_CD./iteration_array_CD,'-^')

xlabel('Number of states')
ylabel('Time/iteration [s]')
legend('sequential','sequential (bad initial guess)','coordinate-descent','Location','northwest')
xticks((2:Nmax)*2) 

%% get number of linear constraints and decision variables (sequential)
n_con = nan(length(solverStats_Seq),1);
n_dec = nan(length(solverStats_Seq),1);
for j = 1:length(solverStats_Seq)
   n_con(j) = solverStats_Seq{j}{1}{end}.conic.size_A.size(1);
   n_dec(j) = solverStats_Seq{j}{1}{end}.conic.size_A.size(2);
end



