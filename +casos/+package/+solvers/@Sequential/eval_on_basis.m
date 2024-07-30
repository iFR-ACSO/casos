function argout = eval_on_basis(obj,argin)

    % Evaluate sequential SOS problem.
    import casos.package.UnifiedReturnStatus
    
    args = argin;
      
    % initialize iteration info
    info = cell(1,obj.opts.max_iter);
    
    % start sequential SOS
    printf(obj.log,'debug','%-15s%-15s%-20s%-20s%-20s%-15s%-10s\n', 'Iteration', 'cost', '||pr||_inf', '||du||_inf', 'squared_L2_vio','alpha','||dLdx||');
    printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');

    % first initial guess must be provided by user! 
    xi_k = args{1};

    % parameter from nonlinear problem
    p0   = args{2};
    
    % initialize Hessian/BFGS approximation
    Bk = eye(obj.sizeHessian(1));

    dual_k = [];
    
    % initilize filter
	if ~isempty(obj.projConPara)

		solPara_proj   = obj.projConPara('p',[obj.xk1fun(xi_k,p0);p0]);
		curr_conVio    = full(solPara_proj.f);

		curr_cost   = inf; %full(casadi.DM(full(obj.f(xi_k,p0))));
       

    else

        curr_conVio = 0;
        curr_cost   = inf;

	end
    
    counter = 0;
    
    feasResdone = 0;

	obj.Filter =  obj.Filter.initializeFilter(curr_conVio ,curr_cost);

    % measure overall solve time
	measTime_seqSOS_in = tic; 

	% solve nonconvex SOS problem via sequence of convex SOS problem
    for i = 0:obj.opts.max_iter
        
        % measure time per iteration
        measTime_seqSOS__iter_in = tic;

        % initial guesss; actually not used but just for completness
        args{1} = xi_k;
    
        % set parameter to convex problem
        args{2}  = [p0; xi_k; Bk(:)];

        % adjust bounds
        args{3}  = argin{3} - xi_k;
        args{4}  = argin{4} - xi_k;

        % evaluate convex SOS problem
        sol = eval_on_basis(obj.sossolver, args);
           
        % store iteration info
        info{i+1}.subProb_stats = obj.sossolver.get_stats;
        
        % check if subproblem is feasible
        switch ( obj.sossolver.get_stats.UNIFIED_RETURN_STATUS)
            case UnifiedReturnStatus.SOLVER_RET_SUCCESS     % optimal soultion found

                    d_star    = sol{1};
                    dual_star = sol{5};

            case UnifiedReturnStatus.SOLVER_RET_UNKNOWN % problem status unknown or ill-posed

                    d_star    = sol{1};
                    dual_star = sol{5};

            case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE
                   % go to feasibility restoration phase

                   % solve feasibility problem i.e. find a feasible
                   % solution that just lies on the constraint manifold



					printf(obj.log,'debug',['Subproblem infeasible in iteration ' num2str(i)   '. Go to feasibility restoration.\n']);
    
					% prepare polynomial input to high-level solver
					polySol = obj.xk1fun(xi_k,p0);

					sol_feas_res = obj.solver_feas_res('x0',[casos.PS(polySol) ; obj.s0 ],...
													   'p', [casos.PS(polySol);0;inf;inf]); 


					xi_feas      = poly2basis(sol_feas_res.x(1:obj.size_x));

					[obj.Filter,~,accept_FeasResStep] = obj.Filter.updateFilter(1,...
																				full(casadi.DM(xi_feas)), ...
																				curr_cost, ....
																				curr_conVio,....
																				full(casadi.DM(full(obj.f(xi_feas,p0)))),...
																				full(sol_feas_res.f),....
																				full(casadi.DM(full(obj.nabla_xi_f(xi_feas,p0))))');

            if accept_FeasResStep
               % solve new iterate

                % just fill info struct for consitency
                info{i+1}.filter_stats.measTime_proj_out      = 0;
                info{i+1}.filter_stats.alpha_k                = 1;
                info{i+1}.filter_stats.measTime               = 0;
                info{i+1}.seqSOS_common_stats.delta_prim      = nan;
                info{i+1}.seqSOS_common_stats.delta_dual      = nan;
                info{i+1}.seqSOS_common_stats.conViol         = curr_conVio ;
                info{i+1}.seqSOS_common_stats.gradLang        = nan;

                info{i+1}.seqSOS_common_stats.solve_time_iter = toc(measTime_seqSOS__iter_in);

				xi_k  = xi_feas;
				printf(obj.log,'debug', 'Feasibility restoration iterate accepted to filter. Continue original problem \n');
				continue
                
				
            else
			
                info{i+1}.seqSOS_common_stats.solve_time = toc(measTime_seqSOS_in);
                error('Problem seems locally infeasible!')
				
            end
             

            otherwise   % problem status unknown or ill-posed

                d_star    = sol{1};
                dual_star = sol{5};
        end

       

        %% backtracking line search filter 
        alpha_max            = obj.opts.alpha_max;
   
        tau                  = obj.opts.tau;
        
        alpha_k = alpha_max;
		
        measTime_filter_in = tic;
        while true

            % compute alpha min 
            dirVecDeriv = full(casadi.DM(full(obj.nabla_xi_f(xi_k,p0)))*d_star);
            if dirVecDeriv < 0
                s_phi = 2.3;
                s_theta =  1.1;
                gamma_phi = 1;
                delta = 1;
                alpha_min = 0.05*min([ 1e-5, ( gamma_phi*curr_conVio )/abs(dirVecDeriv ), (delta*curr_conVio^s_theta )/(abs(dirVecDeriv)^s_phi) ]);

            else

                alpha_min = 1e-5; % same value as gamma_theta in filter

            end

            if alpha_k < alpha_min

					printf(obj.log,'debug',['Step length below minimum step length in iteration ' num2str(i)   '. Go to feasibility restoration.\n']);
    
					% prepare polynomial input to high-level solver
					polySol = obj.xk1fun(xi_k,p0);
                    measTime_feasRes_in = tic;

                    eps_opt = max ([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0)))),(new_conVio)])/obj.opts.optTol;

                    if eps_opt > 1
                        eps_opt = 1;
                    end

                    obj.Filter.augmentFilter(new_cost,new_conVio);

                    zeta_val = 1-eps_opt^2;
                        

					sol_feas_res = obj.solver_feas_res('x0',[casos.PS(polySol) ; obj.s0 ],...
													   'p', [casos.PS(polySol);zeta_val;inf;curr_conVio]); 


					xi_feas      = poly2basis(sol_feas_res.x(1:obj.size_x));

                   measTime_Proj_in = tic;
                    solPara_proj = obj.projConPara('p',[obj.xk1fun(xi_feas,p0);p0]);
                    new_conVio   = full(solPara_proj.f);
                   info{i+1}.filter_stats.measTime_proj_out = toc(measTime_Proj_in);

					[obj.Filter,~,accept_FeasResStep] = obj.Filter.updateFilter(1,...
																				full(casadi.DM(xi_feas)), ...
																				curr_cost, ....
																				curr_conVio,....
																				full(casadi.DM(full(obj.f(xi_feas,p0)))),...
																				new_conVio,....
																				full(casadi.DM(full(obj.nabla_xi_f(xi_feas,p0))))');

                    measTime_feasRes_out = toc(measTime_feasRes_in);

            if accept_FeasResStep
               % solve new iterate

                % just fill info struct for consitency
                % info{i+1}.filter_stats.measTime_proj_out      = 0;
                info{i+1}.filter_stats.alpha_k                = alpha_k;
                info{i+1}.filter_stats.measTime               = measTime_feasRes_out;
                info{i+1}.seqSOS_common_stats.delta_prim      = nan;
                info{i+1}.seqSOS_common_stats.delta_dual      = nan;
                info{i+1}.seqSOS_common_stats.conViol         = new_conVio ;
                info{i+1}.seqSOS_common_stats.gradLang        = full(casadi.DM(full(obj.nabla_xi_L_norm(xi_feas,dual_star,p0))));

                info{i+1}.seqSOS_common_stats.solve_time_iter = toc(measTime_seqSOS__iter_in);
      
				xi_k  = xi_feas;

				printf(obj.log,'debug', 'Feasibility restoration iterate accepted to filter. Continue original problem \n');
                feasResdone = 1;
				break
                
				
            else
			
                info{i+1}.seqSOS_common_stats.solve_time = toc(measTime_seqSOS_in);
                error('Problem seems locally infeasible!')
				
            end
            end
                
            % new-iterate
            xi_k1 = xi_k + alpha_k*d_star;

            % check cost, gradient of cost and constraint violation
            new_cost     = full(casadi.DM(full(obj.f(xi_k1,p0))));
            nabla_xi_f   = full(casadi.DM(full(obj.nabla_xi_f(xi_k1,p0))))';
            
            % check constraint violation (projection onto SOS cone)
            if ~isempty(obj.projConPara)
                
                % measure time to compute projection/current constraint
                % violation
                measTime_Proj_in = tic;
                    solPara_proj = obj.projConPara('p',[obj.xk1fun(xi_k1,p0);p0]);
                    new_conVio   = full(solPara_proj.f);
                info{i+1}.filter_stats.measTime_proj_out = toc(measTime_Proj_in);

            else

                new_conVio = 0;

            end

                   
            % check filter acceptance
            [obj.Filter,goto_SOC,accept_new_iter] = obj.Filter.updateFilter(alpha_k,...
                                                                            full(d_star), ...
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
                       args_solver{2} = [p0; xi_k; d_corr_full; Bk(:)];
        
                        % evaluate subproblem of SOC
                       sol_soc = eval_on_basis(obj.solver_soc, args_solver);
                       
                       % extract adjusted search direction
                       d_corr = sol_soc{1};
        
                       % compute new corrected search direction
                       d_corr_full = d_corr_full + d_corr;
                       
                       % compute new trial point based on corrected
                       % direction; alpha = 1;
                       x_soc_trial = xi_k + alpha_k*d_corr_full;  
        
                       % check cost, gradient of cost 
                       new_cost     = full(casadi.DM(full(obj.f(x_soc_trial,p0))));
                       nabla_xi_f   = full(casadi.DM(full(obj.nabla_xi_f(x_soc_trial,p0))))';
                       
                       % check constraint violation (projection onto SOS cone)
                       solPara_proj = obj.projConPara('p',[obj.xk1fun(x_soc_trial,p0);p0]);
                       new_conVio   = full(solPara_proj.f);
                           
                       % check if new trial point can be accepted to filter
                       % and fulfills sufficient decrease conditions
                       [obj.Filter,repeat,accept_new_iter_soc] = obj.Filter.updateFilter(alpha_k,...
                                                                                        full(d_star), ...
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
                                % counter
                                 iter_soc  = iter_soc  + 1;
        
                                curr_conVio = new_conVio;
                           end
                           
                    end % --- end of while SOC ---
                    info{i+1}.filter_stats.measTime_SOC_out = toc(measTime_SOC_in);
        
            else
                    % not acceptable to filter, adjust step length
                    alpha_k = alpha_k*tau ;
                
            end
        
      

        end % --- end while ---

        if feasResdone 
            feasResdone = 0;
            continue
        end

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
			delta_dual_double = norm(full( (dual_k1 - dual_k)),inf)/norm(full(dual_k1),inf);
        
        
			printf(obj.log,'debug','%-15d%-15e%-20e%-20e%-20e%-15.4f%-18e\n',...
                 i, new_cost , delta_xi_double, delta_dual_double, (new_conVio) , alpha_k , (full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_k1,p0)))))    );



			% printf(obj.log,'debug','%-15d%-15e%-20e%-20e%-20e%-15.4f%-18e\n',...
            %      i, new_cost , delta_xi_double, delta_dual_double, sqrt(new_conVio) , alpha_k , norm(full(casadi.DM(full(obj.nabla_xi_L(xi_k1,dual_k1,p0)))),inf)/norm(full(dual_k1),inf)    );

      

            % store common optimization data
            info{i+1}.seqSOS_common_stats.delta_prim  = delta_xi_double;
            info{i+1}.seqSOS_common_stats.delta_dual  = delta_dual_double;
            info{i+1}.seqSOS_common_stats.conViol     = new_conVio;
            info{i+1}.seqSOS_common_stats.gradLang    = full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_star,p0))));


        
        else
        
            delta_xi_double   = norm( full( casadi.DM(xi_k)    ), inf ); 
            delta_dual_double = norm( full( dual_star ), inf );
            
			printf(obj.log,'debug','%-15d%-15e%-20e%-20e%-20e%-15.4f%-18e\n',...
                 i, new_cost, delta_xi_double, delta_dual_double, (new_conVio), alpha_k , (full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0)))))   );
            
            % store common optimization data
            info{i+1}.seqSOS_common_stats.delta_prim  = delta_xi_double;
            info{i+1}.seqSOS_common_stats.delta_dual  = delta_dual_double;
            info{i+1}.seqSOS_common_stats.conViol     = new_conVio;
            info{i+1}.seqSOS_common_stats.gradLang    = full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0))));
           
        
        end % --- end of display output ---
            
        %% check convergence criteria
        optimality_flag =     max ([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0)))),(new_conVio)]) <= obj.opts.optTol;
        

        % (full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_k1,p0)))))/obj.opts.optTol
        % check if solution stays below tolerance for a certain number of iterations --> solved to acceptable level
        if max ([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0)))),(new_conVio)]) <= obj.opts.accTol

            counter = counter + 1;

        else

            counter = 0;

        end % --- end of counter ---

        
        %% termination criteria and output preparation
        if i ~= obj.opts.max_iter && optimality_flag
                 
                % store iteration info
                info{i+1}.seqSOS_common_stats.solve_time_iter  = toc(measTime_seqSOS__iter_in);
                info{i+1}.seqSOS_common_stats.solve_time = toc(measTime_seqSOS_in);
                info(i+2:end) = [];
                obj.info.iter = info;
               
                % adjust last solution similar to iteration i.e. overwrite optimization solution 
                sol{1} = xi_k1;
                sol{5} = dual_k1;
               
                argout = sol;
        
                printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
                printf(obj.log,'debug','Solution status: Optimal solution found\n'); 
                printf(obj.log,'debug',['Solution time: ' num2str(toc(measTime_seqSOS_in)) ' s\n']); 
                
                % terminate
                return

        elseif i ~= obj.opts.max_iter && counter == obj.opts.noAccIter
    
                % store iteration info
                info{i+1}.seqSOS_common_stats.solve_time_iter  = toc(measTime_seqSOS__iter_in);
                info{i+1}.seqSOS_common_stats.solve_time       = toc(measTime_seqSOS_in);
				info(i+2:end) = [];
				obj.info.iter = info;
               
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
        s = xi_k1-xi_k; 
        y = casadi.DM( full(obj.nabla_xi_L(xi_k1,dual_k1,p0) - obj.nabla_xi_L(xi_k,dual_k1,p0)))';
        
        % prepared hessian approximation for next iteration
        Bk = dampedBFGS(obj,Bk,s,y);
        

        
        %% new solution becomes new linearization point
        xi_k   = xi_k1;
        dual_k = dual_k1;
        
        info{i+1}.seqSOS_common_stats.solve_time_iter  = toc(measTime_seqSOS__iter_in);
    end % --- end for-loop ---

    argout = sol;

end % --- end of function ---
