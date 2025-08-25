function argout = eval_on_basis(obj,argin)
% Solve sequential SOS problem.
measTime_seqSOS_in = tic;
import casos.package.UnifiedReturnStatus

%% print output with current problem and setting
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug',['\t\t Ca' char(931) 'oS-Nonlinear sum-of-squares optimization suite \n']);
printf(obj.log,'debug','\t\t Sequential Quadratic Sum-of-Squares Solver v1.0, April''25 \n');
printf(obj.log,'debug','\t\t GNU GENERAL PUBLIC LICENSE \n');
printf(obj.log,'debug','\t\t Institute of Flight Mechanics and Control, Univ. of Stuttgart \n');
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug','Problem:\n Decision variables = %d (coefficents) \n SOS constraints \t= %d \n Linear constraints = %d\n',obj.display_para.no_decVar   , obj.display_para.no_sosCon, 0 );
printf(obj.log,'debug','Settings:\n Solver:\t%s\n max_iter:\t%d \n Con. Violation Check: \t%s\n', obj.display_para.solver,obj.opts.max_iter,'signed-distance');
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug','%-8s%-15s%-15s%-15s%-15s%-10s%-10s\n', 'iter', 'obj', ['||' char(916) 'x'  '||_inf'], ['||' char(916) char(955)  '||_inf'], [char(920) '(x_k)'],char(945),'||dLdx||');
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');

initTimeMeas = tic;
% input arguments
args = argin;

% initial guess (user)
x_k = args{1};

% parameter from nonlinear problem
p0   = args{2};

% initialize iteration info struct
info.single_iterations = cell(1,obj.opts.max_iter);

% initialze dual variables
dual_k = zeros(obj.init_para.no_dual_var,1);

if strcmp(obj.opts.Hessian_init,'Identity')

    % initialize BFGS-matrix = identity x scaling
    Bk =  eye(obj.init_para.size_B)*obj.opts.scale_BFGS0;

elseif strcmp(obj.opts.Hessian_init,'Analytical')

    % initialize Hessian approximation with regularization
    H = full(obj.eval_Hessian(x_k,p0,dual_k));

    Bk = casos.package.solvers.SequentialCommon.regularize_Hessian(H);

end

if strcmp(obj.opts.conVioCheck,'signed-distance')
    % initialize filter and first iterate
    args_conVio     =  args;
    args_conVio{2}  =  [p0; x_k];
    args_conVio{3}  = -inf(obj.init_para.conVio.no_con,1);
    args_conVio{4}  =  inf(obj.init_para.conVio.no_con,1);

    % compute constraint violation of initial guess
    sol_convio = eval_on_basis(obj.solver_conVio, args_conVio);

    theta_x_k = full(max(0,max(sol_convio{1})));
else
    g_val = full(obj.eval_constraintSamples(x_k,obj.opts.userSample));

    min_vio = min(min(g_val));
    if min_vio < 0
        theta_x_k = abs(min_vio);
    else
        theta_x_k = 0;
    end
end
% evalaute cost function of nonlinear problem
f_x_k     = full(obj.eval_cost(x_k,p0));

% add some boundary to it
filter = [max(1,theta_x_k)*10, inf];

iter = 1;

% just for display output in first iteration
delta_xi_double   = norm(full( casadi.DM(x_k)),inf);
delta_dual_double = norm(full( (dual_k)),inf);
alpha_k           = 1;

counterAcceplvl  = 0;
feasibility_flag = 1;

info.initTime = toc(initTimeMeas);

while iter <= obj.opts.max_iter

    % display output header every 10th iteration
    if ~mod(iter,10) && iter > 0
        printf(obj.log,'debug','%-8s%-15s%-15s%-15s%-15s%-10s%-10s\n', 'iter', 'obj', ['||' char(916) 'x'  '||_inf'], ['||' char(916) char(955)  '||_inf'], [char(920) '(x_k)'],char(945),'||dLdx||');
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
    end

    % display output current iterate
    printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
        iter-1, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , full(casadi.DM( full(obj.eval_gradLag(x_k,p0,dual_k)) ))   );


    %% check convergence (first-order optimality)
    if full(obj.eval_gradLag(x_k,p0,dual_k))  <= full((obj.opts.tolerance_opt*max(1,abs(f_x_k)) + obj.eval_gradLag2(x_k,p0,dual_k))/delta_xi_double)  ... % scaled optimality
            && theta_x_k <= obj.opts.tolerance_con && feasibility_flag == 0                                                                                     % constraint violation

        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Optimal solution found\n');
        solveTime = toc(measTime_seqSOS_in);
        printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);
        feasibility_flag = 1;
        info.solutionStatus = 'Optimal solution';
        break

    end

    % check if current iterate is solved to an acceptable level
    if full(obj.eval_gradLag(x_k,p0,dual_k))  <= full((10*obj.opts.tolerance_opt*max(1,abs(f_x_k)) + obj.eval_gradLag2(x_k,p0,dual_k))/delta_xi_double)  ... % scaled optimality ...
            && theta_x_k <= obj.opts.tolerance_con

        % increase counter
        counterAcceplvl = counterAcceplvl +1;

        % we have 5 in a row, seems we solved it to an acceptable level/almost optimal
        if counterAcceplvl >= obj.opts.almostOptCount
            printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
            printf(obj.log,'debug','Solution status: Almost optimal.\n');
            solveTime = toc(measTime_seqSOS_in);
            printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
            printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
            printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);

            feasibility_flag = 1;

            info.solutionStatus = 'Alomst Optimal solution.';
            break
        end
    else
        % if current iterate is not solved to an acceptable level, reset counter
        counterAcceplvl = 0;
    end

    feasibility_flag = 0;

    % compute solution of current iterate: solve Q-SDP, perform linesearch and update BFGS
    [sol_iter,sol_qp,feas_res_flag,info.single_iterations,obj,filter,Bk] = do_single_iteration(obj, ...
        iter,...
        x_k,...
        dual_k,...
        theta_x_k,...
        f_x_k,...
        Bk,...
        p0,...
        args, ...
        filter, ...
        info.single_iterations);

    info.single_iterations{iter}.timeStats.feasResTime = 0;
    % Invoke feasibility restoration if necessary
    if feas_res_flag >= 1 && ...              % restoration requested
            theta_x_k >  obj.opts.tolerance_con    % is the constraint violation larger then tolerance; otherwise restoration does not make sense


        % parameter (just for readabilty)
        gamma_theta = obj.opts.filter_struct.gamma_theta ;
        gamma_phi   = obj.opts.filter_struct.gamma_phi ;

        % augment filter with iterate where restoration is invoked; this point
        % shall be avoided in the future
        filter = [filter;[theta_x_k*(1-gamma_theta), f_x_k - gamma_phi*theta_x_k]];

        tmptimeStats = info.single_iterations{iter}.timeStats;
        measFeasResTime = tic;
        % feasibility restoration phase
        [sol_iter,sol_qp_feas,filter,feasibility_flag,info.single_iterations] = obj.feas_res_solver.eval_extended(filter,x_k,theta_x_k,p0,obj,info.single_iterations,iter);

        info.single_iterations{iter}.timeStats = tmptimeStats;
        info.single_iterations{iter}.timeStats.feasResTime = toc(measFeasResTime);

        if feasibility_flag >=1

            % continue from restored solution
            delta_xi_double   = norm(full( sol_iter.x_k1 - x_k) ,inf);
            delta_dual_double = norm(full( dual_k ),inf);

            x_k                 = sol_iter.x_k1;
            sol_iter.dual_k1    = dual_k;            % we keep dual soultion from previous iterate
            theta_x_k           = sol_iter.theta_x_k1;
            f_x_k               = sol_iter.f_x_k1;
            alpha_k             = sol_iter.alpha_k;

        else
            % could not recover return iterate where restoration was
            % invoked; return point where invoked
            sol    = sol_qp;

            sol{1} = x_k;
            sol{5} = dual_k;

            break
        end

    elseif feas_res_flag >= 1 &&  theta_x_k < obj.opts.tolerance_con
        % we are already better then treshhold i.e. restoration can't
        % make it better
        sol    = sol_qp;

        % if one, means Q-SDP is infeasible, thus sol.iter should be
        % zero hence return "old" iterate
        if feas_res_flag == 1
            sol{1} = x_k;
            sol{5} = dual_k;
        else
            sol{1} = sol_iter.x_k1;
            sol{5} = sol_iter.dual_k1;
        end

        % just a final check
        if full(obj.eval_gradLag(x_k,p0,dual_k))  <= obj.opts.tolerance_opt*max(1,full(obj.eval_gradLag(x_k,p0,dual_k))) ...
                && theta_x_k <= obj.opts.tolerance_con

            printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
            printf(obj.log,'debug','Solution status: Optimal solution found\n');
            solveTime = toc(measTime_seqSOS_in);
            printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
            printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
            printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);
            feasibility_flag = 1;

        else
            % stalled
            feasibility_flag = -2;

        end

        break



    else % continue with normal iteration

        % for display output
        delta_xi_double   = norm(full( sol_iter.x_k1 - x_k) ,inf);
        delta_dual_double = norm(full( sol_iter.dual_k1 - dual_k ),inf);

        % set values for next iterate
        x_k       = sol_iter.x_k1;
        dual_k    = sol_iter.dual_k1;
        theta_x_k = sol_iter.theta_x_k1;
        f_x_k     = sol_iter.f_x_k1;
        alpha_k   = sol_iter.alpha_k;

    end

    % return last solution
    if ~isempty(sol_qp)
        sol    = sol_qp;
    elseif ~isempty(sol_qp_feas)
        sol = sol_qp_feas;
    else
        sol = [];
    end

    % exchange the primal/dual solution of QP with the accepted iterate
    sol{1} = sol_iter.x_k1;
    sol{2} = sol_iter.f_x_k1;

    % dual variables estimated by underlying quadratic SDP
    sol{4} = sol_qp{4};
    sol{5} = sol_iter.dual_k1;

    iter = iter + 1;

end % end of while loop

% final display output for user
if iter >= obj.opts.max_iter

    % display output header
    if ~mod(iter,10) && iter > 0
        printf(obj.log,'debug','%-8s%-15s%-15s%-15s%-15s%-10s%-10s\n', 'iter', 'obj', ['||' char(916) 'x'  '||_inf'], ['||' char(916) char(955)  '||_inf'], [char(920) '(x_k)'],char(945),'||dLdx||');
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
    end

    % display output current iterate
    printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
        iter-1, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , full(casadi.DM( full(obj.eval_gradLag(x_k,p0,dual_k)) ))   );


    if full(obj.eval_gradLag(x_k,p0,dual_k))  <= obj.opts.tolerance_opt*max(1,full(obj.eval_gradLag(x_k,p0,dual_k))) ...
            && theta_x_k <= obj.opts.tolerance_con && feasibility_flag == 0

        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Optimal solution found\n');
        solveTime = toc(measTime_seqSOS_in);
        printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);

        info.solutionStatus = 'Optimal Solution';

    else

        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Maximum number of iterations reached\n');
        solveTime = toc(measTime_seqSOS_in);
        printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);

        info.solutionStatus = 'Max. Iterations';
    end
end

if feasibility_flag == 0 && iter < obj.opts.max_iter

    % print last iterate
    printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
        iter, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , full(casadi.DM( full(obj.eval_gradLag(x_k,p0,dual_k)) ) )   );

    printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
    printf(obj.log,'debug','Solution status: Problem infeasible.\n');
    solveTime = toc(measTime_seqSOS_in);
    printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
    printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
    printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);

    info.solutionStatus = 'Problem infeasible';

elseif feasibility_flag == -1 && iter < obj.opts.max_iter


    % print last iterate
    printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
        iter, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , full(casadi.DM( full(obj.eval_gradLag(x_k,p0,dual_k)) ) )   );


    printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
    printf(obj.log,'debug','Solution status: Solver stalled in feasibility restoration.\n');
    solveTime = toc(measTime_seqSOS_in);
    printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
    printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
    printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);

    info.solutionStatus = 'Solver stalled';

elseif feasibility_flag == -2 && iter < obj.opts.max_iter

    printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
    printf(obj.log,'debug','Solution status: Feasible solution. Restoration not invoked because already feasible.\n');
    solveTime = toc(measTime_seqSOS_in);
    printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
    printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
    printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);


    info.solutionStatus = 'Feasible solution';
end

info.single_iterations = info.single_iterations(1:iter-1);

info.totalSolveTime    = solveTime;
info.solverBuildTime   = obj.display_para.solver_build_time;
info.iterations        = iter-1; % we do not count the initial guess
info.timePerIteration  = info.totalSolveTime/info.iterations;

% assign output
argout = sol;

% store iteration info
obj.info = info;

end % end of function
