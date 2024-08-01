function [] = plotSolverStats(stats_struct)

iterNo = length(stats_struct.iter);

%% extract data
% initializes array
alpha               = zeros(iterNo,1);
lineSearchTime      = zeros(iterNo,1);
proj_lineSearch     = zeros(iterNo,1);
solve_time_iter     = zeros(iterNo,1);
delta_dual          = zeros(iterNo,1);
delta_primal        = zeros(iterNo,1);
conVio              = zeros(iterNo,1);
grad_Langrangian    = zeros(iterNo,1);

% extract data from every single iteration
for k = 1:iterNo
        
    % get current iteration
    currIter_struct = stats_struct.iter{k};

    % extract filter data 
    alpha(k)            = currIter_struct.filter_stats.alpha_k;
    lineSearchTime(k)   = currIter_struct.filter_stats.measTime;
    proj_lineSearch(k)  = currIter_struct.filter_stats.measTime_proj_out ;

    % extract common optimization data 
    delta_dual(k)       = currIter_struct.seqSOS_common_stats.delta_dual;
    delta_primal(k)     = currIter_struct.seqSOS_common_stats.delta_prim;
    conVio(k)           = currIter_struct.seqSOS_common_stats.conViol;
    grad_Langrangian(k) = currIter_struct.seqSOS_common_stats.gradLang;
    solve_time_iter(k)  = currIter_struct.seqSOS_common_stats.solve_time_iter;

    if k == iterNo
        solve_time  = currIter_struct.seqSOS_common_stats.solve_time;
    end
    xTickLabels{k} = num2str(k);
end



%% plotting
figure('Name','Plot Stats Filter Linesearch')
set(gcf,'Position',[573,329.6666666666666,986,398])
subplot(1,3,1)
plot(0:iterNo-1,alpha,'*-')
axis([0 iterNo-1 0 1.1])
xlabel('Iteration')
ylabel('step length \alpha')

subplot(1,3,2)
semilogy(0:iterNo-1,lineSearchTime)
axis([0 iterNo-1 1/100*min(lineSearchTime) 100*max(lineSearchTime)])
xlabel('Iteration')
ylabel('Meas. Time Linesearch in [s]')

subplot(1,3,3)
semilogy(0:iterNo-1,proj_lineSearch)
axis([0 iterNo-1 1/100*min(proj_lineSearch) 100*max(proj_lineSearch)])
xlabel('Iteration')
ylabel('Meas. Time constraint violation in [s]')


figure('Name','Plot Common Optimization Related Data')
set(gcf,'Position',[573,329.6666666666666,986,398])
subplot(1,3,1)
plot(0:iterNo-1,conVio,'*-')
xlabel('Iteration')
ylabel('Constraint Violation \theta(\xi,k)')

subplot(1,3,2)
plot(0:iterNo-1,grad_Langrangian,'*-')
xlabel('Iteration')
ylabel('||\nabla_{\xi,k} L||')

subplot(1,3,3)
semilogy(0:iterNo-1,solve_time_iter)
axis([0 iterNo-1 1/100*min(solve_time_iter) 100*max(solve_time_iter)])
xlabel('Iteration')
ylabel('Meas. Time per Iteration in [s]')



end