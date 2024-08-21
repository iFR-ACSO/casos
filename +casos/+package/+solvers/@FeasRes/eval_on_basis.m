function argout = eval_on_basis(obj,argin)

    % Evaluate sequential SOS problem.
    import casos.package.UnifiedReturnStatus
    
    args = argin;
      
    % initialize iteration info
    info = cell(1,obj.opts.max_iter);
    
    % start sequential SOS
    printf(obj.log,'debug','%-15s%-15s%-20s%-20s%-20s%-15s%-10s\n', 'Iteration', 'cost', '||pr||_inf', '||du||_inf', '||conVio||_2','alpha','||dLdx||');
    printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');

    % first initial guess must be provided by user! 
    xi_k = args{1};

    % parameter from nonlinear problem
    p0   = args{2};
    
    % initialize Hessian/BFGS approximation
    Bk = eye(obj.sizeHessian(1));

    dual_k = [];
    
    counter_acceptLvl = 0;
   
    % initilize filter with constrain violation of last iterate
    curr_conVio = full(casadi.DM(p0(2)));
    curr_cost   = inf;

	obj.Filter =  obj.Filter.initializeFilter( curr_conVio ,curr_cost);

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

            otherwise %case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE
                   

                % return starting point
                xi_feas   = xi_k;
                % dual_star = sol{5};

                info{i+1}.seqSOS_common_stats.solve_time = toc(measTime_seqSOS_in);
 
                info{i+1}.filter_stats.measTime_proj_out      = 0;
                info{i+1}.filter_stats.alpha_k                = 1;
                info{i+1}.filter_stats.measTime               = nan;
                info{i+1}.seqSOS_common_stats.delta_prim      = nan;
                info{i+1}.seqSOS_common_stats.delta_dual      = nan;
                info{i+1}.seqSOS_common_stats.conViol         = nan ;
                info{i+1}.seqSOS_common_stats.gradLang        = nan;
                info{i+1}.seqSOS_common_stats.solve_time_iter = toc(measTime_seqSOS__iter_in);

                obj.info.UNIFIED_RETURN_STATUS_feas = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;

                info(i+2:end) = [];
                obj.info.iter = info;
               
                % adjust last solution similar to iteration i.e. overwrite optimization solution 
                sol{1} = xi_feas;
                % sol{5} = dual_k1;
               
                argout = sol;
        
                printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
                printf(obj.log,'debug','Solution status: Solver infeasible in restoration phase!\n'); 
                printf(obj.log,'debug',['Solution time: ' num2str(toc(measTime_seqSOS_in)) ' s\n']); 
                
                % terminate
                return
           
				          
        end

       

        %% backtracking line search filter 

        % extract parameter (just for better code reading)
        alpha_max            = obj.opts.alpha_max;
        tau                  = obj.opts.tau;
        alpha_k              = alpha_max;
		
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
            
  
                % measure time to compute projection/current constraint
                % violation
                new_conVio = sqrt(full(casadi.DM(obj.proj_conFun(xi_k1,p0)))) ;
                % measTime_Proj_in = tic;
                %     solPara_proj = obj.projConPara('p',[obj.xk1fun(xi_k1,p0)]);
                %     new_conVio   = sqrt(full(solPara_proj.f));
                % info{i+1}.filter_stats.measTime_proj_out = toc(measTime_Proj_in);



                  
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
                       solPara_proj = obj.projConPara('p',obj.xk1fun(x_soc_trial,p0));
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
                                % counter_acceptLvl
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
                    if  max([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))))/obj.opts.optTol,(new_conVio)/obj.opts.conVioTol]) <= 1                                                                                                            
                         break                                                                                                          
                    end
                        
                    % solver stalled
                    if alpha_k < alpha_min
                       
			        
                        % not acceptable to filter after feasibility
                        % restoration --> solver stalled
                        info{i+1}.seqSOS_common_stats.solve_time = toc(measTime_seqSOS_in);
        
                        info{i+1}.filter_stats.measTime_proj_out      = 0;
                        info{i+1}.filter_stats.alpha_k                = alpha_k;
                        info{i+1}.filter_stats.measTime               = measTime_feasRes_out;
                        info{i+1}.seqSOS_common_stats.delta_prim      = nan;
                        info{i+1}.seqSOS_common_stats.delta_dual      = nan;
                        info{i+1}.seqSOS_common_stats.conViol         = new_conVio ;
                        info{i+1}.seqSOS_common_stats.gradLang        = full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))));
                        info{i+1}.seqSOS_common_stats.solve_time_iter  = toc(measTime_seqSOS__iter_in);
                        info(i+2:end) = [];
                        obj.info.iter = info;
                        obj.info.UNIFIED_RETURN_STATUS_feas = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    
                        % adjust last solution similar to iteration i.e. overwrite optimization solution 
                        sol{1} = xi_k1;
                        sol{5} = dual_k1;
                       
                        argout = sol;
                
                        printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
                        printf(obj.log,'debug','Solution status: Feasibility Restoration Stalled\n'); 
                        printf(obj.log,'debug',['Solution time: ' num2str(toc(measTime_seqSOS_in)) ' s\n']); 
                        
                        % terminate
                        return

                    end
            end 
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
        
        
			printf(obj.log,'debug','%-15d%-15e%-20e%-20e%-20e%-15.4f%-18e\n',...
                 i, new_cost , delta_xi_double, delta_dual_double, new_conVio , alpha_k ,...
                 full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))))    );


            % store common optimization data
            info{i+1}.seqSOS_common_stats.delta_prim  = delta_xi_double;
            info{i+1}.seqSOS_common_stats.delta_dual  = delta_dual_double;
            info{i+1}.seqSOS_common_stats.conViol     = new_conVio;
            info{i+1}.seqSOS_common_stats.gradLang    = full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))));
        
        else
            % just for the first iteration
            delta_xi_double   = norm( full( casadi.DM(xi_k)    ), inf ); 
            delta_dual_double = norm( full( dual_star ), inf );
            
			printf(obj.log,'debug','%-15d%-15e%-20e%-20e%-20e%-15.4f%-18e\n',...
                 i, new_cost, delta_xi_double, delta_dual_double, new_conVio, alpha_k ,...
                 full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0))))   );
            
            % store common optimization data
            info{i+1}.seqSOS_common_stats.delta_prim  = delta_xi_double;
            info{i+1}.seqSOS_common_stats.delta_dual  = delta_dual_double;
            info{i+1}.seqSOS_common_stats.conViol     = new_conVio;
            info{i+1}.seqSOS_common_stats.gradLang    = full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0))));
           
        
        end % --- end of display output ---
            
        %% check convergence criteria

        optimality_flag =     max ([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))))/...
                                    obj.opts.optTol,new_conVio/obj.opts.conVioTol]) <= 1;
        

        % check if solution stays below tolerance for a certain number of iterations --> solved to acceptable level
        if max ([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))))/...
                      obj.opts.accTol,new_conVio/obj.opts.conVioTol/10]) <= 1
 
            counter_acceptLvl = counter_acceptLvl + 1;

        else

            counter_acceptLvl = 0;

        end % --- end of counter ---

        
        %% termination criteria and output preparation
        if i ~= obj.opts.max_iter && optimality_flag
                 
                % store iteration info
                info{i+1}.seqSOS_common_stats.solve_time_iter  = toc(measTime_seqSOS__iter_in);
                info{i+1}.seqSOS_common_stats.solve_time       = toc(measTime_seqSOS_in);
                info(i+2:end) = [];
                obj.info.iter = info;
               
                obj.info.UNIFIED_RETURN_STATUS_feas = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
                % adjust last solution similar to iteration i.e. overwrite optimization solution 
                sol{1} = xi_k1;
                sol{5} = dual_k1;
               
                argout = sol;
        
                printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
                printf(obj.log,'debug','Solution status: Optimal solution found\n'); 
                printf(obj.log,'debug',['Solution time: ' num2str(toc(measTime_seqSOS_in)) ' s\n']); 
                
                % terminate
                return

        elseif i ~= obj.opts.max_iter && counter_acceptLvl == obj.opts.noAccIter
    
                % store iteration info
                info{i+1}.seqSOS_common_stats.solve_time_iter  = toc(measTime_seqSOS__iter_in);
                info{i+1}.seqSOS_common_stats.solve_time       = toc(measTime_seqSOS_in);
				info(i+2:end) = [];
				obj.info.iter = info;
               
                obj.info.UNIFIED_RETURN_STATUS_feas = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS; 
                % adjust last solution similar to iteration i.e. overwrite optimization solution 
                sol{1} = xi_k1;
                sol{5} = dual_k1;
               
                argout = sol;
        
                printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
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
        
                printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
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
        
        
        %% new solution becomes new linearization point
        xi_k   = xi_k1;
        dual_k = dual_k1;
        
        % measure time needed for one full iteration
        info{i+1}.seqSOS_common_stats.solve_time_iter  = toc(measTime_seqSOS__iter_in);

    end % --- end for-loop ---

    argout = sol;

end % --- end of function ---
