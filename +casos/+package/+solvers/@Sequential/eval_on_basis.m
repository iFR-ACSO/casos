function argout = eval_on_basis(obj,argin)

% Evaluate sequential SOS problem.
import casos.package.UnifiedReturnStatus

args = argin;

% initialize iteration info
info = cell(1,obj.opts.max_iter);

%% print output with current problem and setting
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug','\t\t Sequential Quadratic Sum-of-Squares Solver v1.0, September''24 \n');
printf(obj.log,'debug','\t\t Institute of Flight Mechanics and Control, Univ. of Stuttgart \n');
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug','Problem:\n Decision variables = %d (coefficents) \n SOS constraints \t= %d \n Linear constraints = %d\n',length(args{1}), 3, 0 );
% printf(obj.log,'debug','Settings:\n Solver:\t%s\n max_iter:\t%d \n Con. Violation Check: \t%s\n', obj.opts.sossol,obj.opts.max_iter,obj.opts.conViolCheck);
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug','%-8s%-15s%-15s%-15s%-15s%-10s%-10s\n', 'iter', 'obj', '||pr||_inf', '||du||_inf', '||conVio||_2','alpha','||dLdx||');
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');

%% prepare problem
% first initial guess must be provided by user!
xi_k = args{1};

% parameter from nonlinear problem
p0   = args{2};

% initialize Hessian/BFGS approximation
Bk = eye(obj.sizeHessian(1));

% To-Do initialize dual variables i.e. default or user value
dual_k = [];

% initalize counter for acceptable level flag
counter_acceptLvl = 0;

% flag for feasibility restoration
feasResdone = 0;

%% initilize filter
if ~isempty(obj.projConPara)
    
    % compute projection of initial guess
    solPara_proj   = obj.projConPara('p',[obj.xk1fun(xi_k,p0);obj.p0poly(p0)]);
    curr_conVio    = full(solPara_proj.f);
    
    curr_cost   = inf;
    
else
    
    curr_conVio = 0;
    curr_cost   = inf;
    
end

% actual initialization
obj.Filter =  obj.Filter.initializeFilter(curr_conVio ,curr_cost);

%% solve problem
% measure overall solve time
measTime_seqSOS_in = tic;

% solve nonconvex SOS problem via sequence of convex SOS problem
for i = 0:obj.opts.max_iter
    
    % measure time per iteration
    measTime_seqSOS__iter_in = tic;
    
    %% prepare input for underlying solver
    % initial guesss; actually not used but just for completness
    args{1} = xi_k;
    
    % set parameter to convex problem
    
    args{2}  = [p0; xi_k; Bk(:);zeros(length(xi_k),1)];
    
    % adjust bounds
    % args{3}  = argin{3} - xi_k;
    % args{4}  = argin{4} - xi_k;
    
    args{3}  = argin{3};
    args{4}  = argin{4};
    %
    %% evaluate convex SOS problem
    sol = eval_on_basis(obj.sossolver, args);
    
    % store iteration info from low level solver
    info{i+1}.subProb_stats = obj.sossolver.get_stats;
    
    %% check if subproblem is feasible
    switch ( obj.sossolver.get_stats.UNIFIED_RETURN_STATUS)
        case UnifiedReturnStatus.SOLVER_RET_SUCCESS     % optimal soultion found
            
            d_star    = sol{1}-xi_k;
            dual_star = sol{5};
            
        case UnifiedReturnStatus.SOLVER_RET_UNKNOWN % feasible solution but not optimal
            
            d_star    = sol{1}-xi_k;
            dual_star = sol{5};
            
        case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE % actually this is otherwise
            
           
                printf(obj.log,'debug',['Subproblem infeasible in iteration ' num2str(i)   '. Feasibility restoration is currently turned off.\n']);
                

                info{i+1}.filter_stats.measTime_proj_out        = 0;
                info{i+1}.filter_stats.alpha_k                  = 1;
                info{i+1}.filter_stats.measTime                 = 0;
                info{i+1}.seqSOS_common_stats.delta_prim        = nan;
                info{i+1}.seqSOS_common_stats.delta_dual        = nan;
                info{i+1}.seqSOS_common_stats.conViol           = curr_conVio ;
                info{i+1}.seqSOS_common_stats.gradLang          = full(casadi.DM(full(obj.nabla_xi_L_norm(xi_feas,dual_star,p0))));
                info{i+1}.seqSOS_common_stats.solve_time_iter   = toc(measTime_seqSOS__iter_in);
                info(i+2:end)                                   = [];
                obj.info.iter                                   = info;
                obj.info.UNIFIED_RETURN_STATUS_seq              = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;

                sol{1} = xi_k;
                sol{5} = dual_k1;
                
                argout = sol;
                
                printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
                printf(obj.log,'debug','Solution status: Solver stalled\n');
                printf(obj.log,'debug',['Solution time: ' num2str(toc(measTime_seqSOS_in)) ' s\n']);
                
                % terminate
                return
                               
    end
    
    
    
    %% backtracking line search filter
    
    % extract parameter (just for code reading)
    alpha_max            = obj.opts.alpha_max;
    tau                  = obj.opts.tau;
    alpha_k              = alpha_max;
    
    % measure time spend in filter
    measTime_filter_in = tic;
    while true
        
        % compute alpha min
        alpha_min = computeMin_stepLength(obj,xi_k,p0,d_star,curr_conVio);
        
        % new-iterate
        xi_k1 = xi_k + alpha_k*d_star;
        
        % check cost, gradient of cost and constraint violation of the
        % new iterate
        new_cost     = full(casadi.DM(full(obj.f(xi_k1,p0))));
        nabla_xi_f   = full(casadi.DM(full(obj.nabla_xi_f(xi_k1,p0))))';
        
        % check constraint violation (projection onto SOS cone)
        if ~isempty(obj.projConPara) && strcmp(obj.opts.conViolCheck,'projection')
            
            % measure time to compute projection/current constraint violation
            measTime_Proj_in = tic;
                % compute exact (in numerical sense) projection/squared distance 
                solPara_proj = obj.projConPara('p',[obj.xk1fun(xi_k1,p0);obj.p0poly(p0)]);

                % L2-norm
                new_conVio   = sqrt( full( solPara_proj.f ) );
            info{i+1}.filter_stats.measTime_proj_out = toc(measTime_Proj_in);
            
        else
            
            new_conVio = 0;
            
            % execute pseudo-projection
            if ~strcmp(obj.opts.conViolCheck,'projection')
                measTime_Proj_in = tic;
                
                % we first evaluate the constraints at the current
                % solution and set it up as a function
                pseudo_proj = to_function(obj.pseudoProj(xi_k1,p0));
                
                % To-Do: sample prep. can be done in build problem
                
                % bring samples into the correct form i.e.
                % pseudo_proj(x_1, x_2, ... , nInputs)
                samples = num2cell(obj.opts.conVioSamp,2);
                
                % evaluate
                new_conVio = min(min(full(pseudo_proj(samples{:}))));
                
                % if min-value is non-negative then there is no
                % violation (based on samples)
                if new_conVio >= 0
                    new_conVio = 0;
                else
                    % just to stick we the definition of filter i.e.
                    % similar to norm only non-negative values
                    new_conVio = abs(new_conVio);
                end
                info{i+1}.filter_stats.measTime_proj_out = toc(measTime_Proj_in);
            end
            
        end
        
        
        % check filter acceptance
        [obj.Filter,goto_SOC,accept_new_iter] = obj.Filter.updateFilter(alpha_k,...
                                                                        full(casadi.DM(d_star)), ...
                                                                        curr_cost, ....
                                                                        curr_conVio,....
                                                                        new_cost,...
                                                                        new_conVio,....
                                                                        nabla_xi_f );
                                                                    
        
        % if acceptable to filter and sufficient decrease are fulfilled
        % accept current iterate
        if accept_new_iter && ~goto_SOC
            break
            
        % second-order correction
        elseif ~accept_new_iter && goto_SOC && alpha_k == alpha_max
            
            % start at current search direction
            d_corr_full = d_star;
            
            iter_soc    = 1;
            
            measTime_SOC_in = tic;
            % perform for maximum number of SOC iterations
            while iter_soc  < obj.opts.soc_max_iter
                
                % adjust parameter for soc
                args_solver    = args;
                
                args_solver{2} = [p0; xi_k; Bk(:);d_corr_full];
                
                % evaluate subproblem of SOC (convex problem)
                sol_soc = eval_on_basis(obj.sossolver, args_solver);
                
                % extract adjusted search direction
                d_corr = sol_soc{1}-xi_k;
                
                % compute new corrected search direction
                d_corr_full = d_corr_full + d_corr;
                
                % compute new trial point based on corrected
                % direction; alpha = 1;
                x_soc_trial = xi_k + alpha_k*d_corr_full;
                
                % check cost, gradient of cost
                new_cost     = full(casadi.DM(full(obj.f(x_soc_trial,p0))));
                nabla_xi_f   = full(casadi.DM(full(obj.nabla_xi_f(x_soc_trial,p0))))';
                
                % check constraint violation (projection onto SOS cone)
                solPara_proj = obj.projConPara('p',[obj.xk1fun(x_soc_trial,p0);obj.p0poly(p0)]);
                new_conVio   = sqrt(full(solPara_proj.f));
                
                % check if new trial point can be accepted to filter
                % and fulfills sufficient decrease conditions
                [obj.Filter,repeat,accept_new_iter_soc] = obj.Filter.updateFilter(alpha_k,...
                                                                                    full(casadi.DM(d_star)), ...
                                                                                    curr_cost, ....
                                                                                    curr_conVio,....
                                                                                    new_cost,...
                                                                                    new_conVio,....
                                                                                    nabla_xi_f );
                
                if accept_new_iter_soc && ~repeat
                    
                    xi_k1  = x_soc_trial;
                    d_star = d_corr_full;
                    break
                    
                elseif ~accept_new_iter_soc && ~repeat
                    % not acceptable to filter, adjust step length
                    alpha_k = alpha_k*tau ;
                    break
                else
                    
                    % if we are here, we start SoC step again
                    % with new corrected search direction
                    
                    % counter while loop
                    iter_soc  = iter_soc  + 1;
                    
                    curr_conVio = new_conVio;
                end
                
            end % --- end of while SOC ---
            info{i+1}.filter_stats.measTime_SOC_out = toc(measTime_SOC_in);
            
        else
            % not acceptable to filter, adjust step length
             alpha_k = alpha_k*tau ;
                
             % means we found an optimal solution, altough it is
             % not acceptable to the filter; maybe a blocking entry
             if  max([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))))/...
                        obj.opts.optTol,(new_conVio)/obj.opts.conVioTol]) <= 1
                  break
             end
            
                % got to feasability restoration if we are below alpha_min
             if alpha_k < alpha_min
                
             
                    % no feasibility restoration phase active
                    printf(obj.log,'debug',['Step length below minimum step length in iteration ' num2str(i)  '. Feasibility restoration is currently turn off.\n']);
                    
                    % prepare output of additional metrics
                    info{i+1}.seqSOS_common_stats.solve_time        = toc(measTime_seqSOS_in);
                    info{i+1}.filter_stats.measTime_proj_out        = 0;
                    info{i+1}.filter_stats.alpha_k                  = alpha_k;
                    info{i+1}.filter_stats.measTime                 = 0;
                    info{i+1}.seqSOS_common_stats.delta_prim        = nan;
                    info{i+1}.seqSOS_common_stats.delta_dual        = nan;
                    info{i+1}.seqSOS_common_stats.conViol           = new_conVio ;
                    info{i+1}.seqSOS_common_stats.gradLang          = full(casadi.DM(full(obj.nabla_xi_L_norm(xi_feas,dual_star,p0))));
                    info{i+1}.seqSOS_common_stats.solve_time_iter   = toc(measTime_seqSOS__iter_in);

                    info(i+2:end) = [];
                    obj.info.iter = info;

                    obj.info.UNIFIED_RETURN_STATUS_seq = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
                    
                    sol{1} = xi_k1;
                    sol{5} = dual_k1;
                    
                    argout = sol;
                    
                    printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
                    printf(obj.log,'debug','Solution status: Solver stalled\n');
                    printf(obj.log,'debug',['Solution time: ' num2str(toc(measTime_seqSOS_in)) ' s\n']);
                    
                    % terminate
                    return
                    
                  
            end % --- end of check for alpha_min ---

        end % --- end of if-else to check current solution ---

    end % --- end while ---

    
    measTime_filter_out = toc(measTime_filter_in);
    
    % store some metrics about filter
    info{i+1}.filter_stats.alpha_k  = alpha_k;
    info{i+1}.filter_stats.measTime = measTime_filter_out;
    
    % prepare new constraint violation/cost for next iteration (filter)
    curr_conVio = new_conVio;
    curr_cost   = new_cost;
    
    % compute dual variables
    if ~isempty(dual_k)
        dual_k1 = dual_k + alpha_k*(dual_star - dual_k);
    else
        dual_k1 = dual_star;
    end
    
   
    %% Prepare display output
    if ~isempty(dual_k)
        
        delta_xi_double   = norm(full( casadi.DM(xi_k1 - xi_k)),inf);
        delta_dual_double = norm(full( (dual_k1 - dual_k)),inf);
        
        
        printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
            i, new_cost , delta_xi_double, delta_dual_double, new_conVio , alpha_k , full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))))   );
        
        
        % store common optimization data
        info{i+1}.seqSOS_common_stats.delta_prim  = delta_xi_double;
        info{i+1}.seqSOS_common_stats.delta_dual  = delta_dual_double;
        info{i+1}.seqSOS_common_stats.conViol     = new_conVio;
        info{i+1}.seqSOS_common_stats.gradLang    = full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))));

    else
        % just for the first iteration
        delta_xi_double   = norm( full( casadi.DM(xi_k)    ), inf );
        delta_dual_double = norm( full( dual_star ), inf );
        
        printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
            i, new_cost, delta_xi_double, delta_dual_double, new_conVio, alpha_k , full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0))))   );
        
        % store common optimization data
        info{i+1}.seqSOS_common_stats.delta_prim  = delta_xi_double;
        info{i+1}.seqSOS_common_stats.delta_dual  = delta_dual_double;
        info{i+1}.seqSOS_common_stats.conViol     = new_conVio;
        info{i+1}.seqSOS_common_stats.gradLang    = full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0))));
        
        
    end % --- end of display output ---
    
    %% check convergence criteria
    
    optimality_flag =     max([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))))/obj.opts.optTol,new_conVio/obj.opts.conVioTol]) <= 1;
    
    
    % check if solution stays below tolerance for a certain number of iterations --> solved to acceptable level
    if max([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))))/obj.opts.accTol,new_conVio/obj.opts.conVioTol/10]) <= 1
        
        counter_acceptLvl = counter_acceptLvl + 1;
        
    else
        
        counter_acceptLvl = 0;
        
    end % --- end of counter ---
    
    %% in case of pseduo-projection check solution flags with exact (in numerical sense) projection to verify
    if (optimality_flag || counter_acceptLvl == obj.opts.noAccIter) && ~strcmp(obj.opts.conViolCheck,'projection')
        if optimality_flag
            
            printf(obj.log,'debug','Solution status: Optimal solution found using pseudo-projection\n');
            printf(obj.log,'debug','Run final check with projection \n');
            
            solPara_proj    = obj.projConPara('p',[obj.xk1fun(xi_k1,p0);p0]);
            new_conVio      = sqrt(full(solPara_proj.f));

            % check flag again
            optimality_flag = max ([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))))/obj.opts.optTol,new_conVio/obj.opts.conVioTol]) <= 1;
            
        else
            printf(obj.log,'debug','Solution status: Solved to acceptable level using pseudo-projection\n');
            printf(obj.log,'debug','Run final check with projection \n');
            
            solPara_proj = obj.projConPara('p',[obj.xk1fun(xi_k1,p0);p0]);
            new_conVio   = sqrt(full(solPara_proj.f));
            
            if max ([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))))/obj.opts.accTol,new_conVio/obj.opts.conVioTol/10]) <= 1
                counter_acceptLvl = counter_acceptLvl + 1;
            end
        end
        
        
    end
    
    
    %% termination criteria and output preparation
    if i ~= obj.opts.max_iter && optimality_flag
        
        % store iteration info
        info{i+1}.seqSOS_common_stats.solve_time_iter  = toc(measTime_seqSOS__iter_in);
        info{i+1}.seqSOS_common_stats.solve_time       = toc(measTime_seqSOS_in);
        
        info(i+2:end) = [];
        obj.info.iter = info;
        
        obj.info.UNIFIED_RETURN_STATUS_seq = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
        
        
        % adjust last solution similar to iteration i.e. overwrite optimization solution
        sol{1} = xi_k1;
        sol{5} = dual_k1;
        
        argout = sol;
        
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Optimal solution found\n');
        printf(obj.log,'debug',['Solution time: ' num2str(toc(measTime_seqSOS_in)) ' s\n']);
        
        % terminate
        return
        
    elseif i ~= obj.opts.max_iter && counter_acceptLvl >= obj.opts.noAccIter
        
        % store iteration info
        info{i+1}.seqSOS_common_stats.solve_time_iter  = toc(measTime_seqSOS__iter_in);
        info{i+1}.seqSOS_common_stats.solve_time       = toc(measTime_seqSOS_in);
        
        info(i+2:end) = [];
        obj.info.iter = info;
        
        obj.info.UNIFIED_RETURN_STATUS_seq = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
        
        % adjust last solution similar to iteration i.e. overwrite optimization solution
        sol{1} = xi_k1;
        sol{5} = dual_k1;
        
        argout = sol;
        
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Solved to acceptable level\n');
        printf(obj.log,'debug',['Solution time: ' num2str(toc(measTime_seqSOS_in)) ' s\n']);
        
        % terminate
        return
        
        
    elseif i == obj.opts.max_iter
        
        % store iteration info
        info{i+1}.seqSOS_common_stats.solve_time_iter  = toc(measTime_seqSOS__iter_in);
        info{i+1}.seqSOS_common_stats.solve_time = toc(measTime_seqSOS_in);
        
        info(i+2:end) = [];
        obj.info.iter = info;
        
        obj.info.UNIFIED_RETURN_STATUS_seq = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
        
        % adjust last solution similar to iteration i.e. overwrite optimization solution
        sol{1} = xi_k1;
        sol{5} = dual_k1;
        
        argout = sol;
        
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Maximum number of iterations reached\n');
        printf(obj.log,'debug',['Solution time: ' num2str(toc(measTime_seqSOS_in)) ' s\n']);
        
        % terminate
        return
        
    end % --- end of output preparation
    
    
    %% Damped BFGS
    
    % update BFGS
    s = full(casadi.DM(xi_k1-xi_k));
    y = full(casadi.DM( full(obj.nabla_xi_L(xi_k1,dual_k1,p0) - obj.nabla_xi_L(xi_k,dual_k1,p0)))');
    
    % prepare hessian approximation for next iteration
    Bk = dampedBFGS(obj,Bk,s,y);
    if obj.opts.debugBFGS
        info{i+1}.BFGS_info.Bk = Bk;
        info{i+1}.BFGS_info.eig = eig(Bk);
        info{i+1}.BFGS_info.con = cond(Bk);
    end

    %% new solution becomes new linearization point
    xi_k   = xi_k1;
    dual_k = dual_k1;
    
    % measure time needed for one full iteration
    info{i+1}.seqSOS_common_stats.solve_time_iter  = toc(measTime_seqSOS__iter_in);
    
end % --- end for-loop ---

argout = sol;

end % --- end of function ---
