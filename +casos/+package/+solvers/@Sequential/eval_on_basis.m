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
    B = eye(obj.sizeHessian(1));
 
    dual_k = [];
    
   if ~isempty(obj.projConPara)

       solPara_proj   = obj.projConPara('p',[obj.xk1fun(xi_k,p0);p0]);
       curr_conVio    = full(solPara_proj.f);

       curr_cost   = full(casadi.DM(full(obj.f(xi_k,p0))));
       

   else
        curr_conVio = 0;
        curr_cost   = inf;
   end
    

   % initialize filter
   obj.Filter =  obj.Filter.initializeFilter(curr_conVio ,curr_cost);

   tic
   % solve nonconvex SOS problem via sequence of convex SOS problem
    for i = 1:obj.opts.max_iter

        % initial guesss; actually not used but just for completness
        args{1} = xi_k;
    
        % set parameter to convex problem
        args{2}  = [p0; xi_k; B(:)];

        % adjust bounds
        args{3}  = argin{3} - xi_k;
        args{4}  = argin{4} - xi_k;

        % evaluate convex SOS problem
        sol = eval_on_basis(obj.sossolver, args);
           
        % store iteration info
        info{i} = obj.sossolver.get_stats;
        
        % check if subproblem is feasible
        switch ( obj.sossolver.get_stats.UNIFIED_RETURN_STATUS)
            case UnifiedReturnStatus.SOLVER_RET_SUCCESS % optimal soultion found

                    d_star    = sol{1};
                    dual_star = sol{5};

            case UnifiedReturnStatus.SOLVER_RET_UNKNOWN % problem status unknown or ill-posed

                    d_star    = sol{1};
                    dual_star = sol{5};

            case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE
                               %% go to feasibility restoration phase
            printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
            printf(obj.log,'debug','Subproblem infeasible: Go to feasibility restoration\n');
    
                % prepare polynomial input to high-level solver
               polySol = obj.xk1fun(xi_k,p0);
    
               sol_feas_res = obj.solver_feas_res('x0',[casos.PS(polySol) ; obj.s0 ],...
                                                  'p',[casos.PS(polySol);0]); 
    
    
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
               xi_k  = xi_feas;
              continue
            else
                error('Problem seems infeasible')
            end
             

            otherwise   % problem status unknown or ill-posed

                    d_star    = sol{1};
                    dual_star = sol{5};
        end

       

        %% backtracking line search filter 
        alpha_max            = 1;
        alpha_min            = 0.0001;
        tau                  = 0.5;
        
        alpha_k = alpha_max;
		
        while alpha_k >= alpha_min
                
            % new-iterate
            xi_k1 = xi_k + alpha_k*d_star;

            % check cost, gradient of cost and constraint violation
            new_cost     = full(casadi.DM(full(obj.f(xi_k1,p0))));
            nabla_xi_f   = full(casadi.DM(full(obj.nabla_xi_f(xi_k1,p0))))';
            
            if ~isempty(obj.projConPara)

                solPara_proj = obj.projConPara('p',[obj.xk1fun(xi_k1,p0);p0]);
                new_conVio   = full(solPara_proj.f);

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

            
            if accept_new_iter && ~goto_SOC
                break
        
            % second-order correction
            elseif ~accept_new_iter && goto_SOC && alpha_k >= alpha_max/2
            
        
                     d_corr_full = d_star;
                     
                     p    = 1;
                     pmax = 5;

                     % perform for maximum number of SOC iterations
                     while p < pmax 
                        
                        % adjust parameter for soc
                        args_solver    = args;
                        args_solver{2} = [p0; xi_k; d_corr_full; B(:)];
        
                         % evaluate subproblem of SOC
                        sol_soc = eval_on_basis(obj.solver_soc, args_solver);
        
                        d_corr = sol_soc{1};
        
                        % compute new corrected search direction
                        d_corr_full = d_corr_full + d_corr;
                        
                        % compute new trial point based on corrected
                        % direction; alpha = 1;
                        x_soc_trial = xi_k + alpha_k*d_corr_full;  
        
                        % check cost, gradient of cost and constraint violation
                        new_cost     = full(casadi.DM(full(obj.f(x_soc_trial,p0))));
                        nabla_xi_f   = full(casadi.DM(full(obj.nabla_xi_f(x_soc_trial,p0))))';
                
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
                                  p = p + 1;
        
                                 curr_conVio = new_conVio;
                            end
                            
                     end % --- end of while SOC ---
        
            else
                    % not acceptable to filter, adjust step length
                    alpha_k = alpha_k*tau ;
                
            end
        
      

        end % --- end while ---

                
            if alpha_k < alpha_min
    
            %% go to feasibility restoration phase
            printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
            printf(obj.log,'debug','step length below minimum step length: Go to feasibility restoration\n');
    
                % prepare polynomial input to high-level solver
               polySol = obj.xk1fun(xi_k1,p0);
                xmax = 1;
               zeta_val     = 1/xmax^2* full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_k1,p0))))^2; % To do: add poly. interpolation
               if zeta_val >= xmax
                   zeta_val = xmax;
               end
               zeta_val = 0.5;
               sol_feas_res = obj.solver_feas_res('x0',[casos.PS(polySol) ; obj.s0 ],...
                                                  'p',  [casos.PS(polySol);zeta_val]); 
    
    
              xi_feas      = poly2basis(sol_feas_res.x(1:obj.size_x));


              
               solPara_proj = obj.projConPara('p',[obj.xk1fun(xi_feas,p0);p0]);
               new_conVio   = full(solPara_proj.f);

              [obj.Filter,~,accept_FeasResStep] = obj.Filter.updateFilter(alpha_k,...
                                                                          full(casadi.DM(xi_feas)), ...
                                                                          curr_cost, ....
                                                                          curr_conVio,....
                                                                          full(casadi.DM(full(obj.f(xi_feas,p0)))),...
                                                                          new_conVio,....
                                                                          full(casadi.DM(full(obj.nabla_xi_f(xi_feas,p0))))');

            if accept_FeasResStep
               % solve new iterate
               xi_k  = xi_feas;
              continue
            else
                error('Problem seems infeasible')
            end
             
    
            end


        % set new iterate for next iteration
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
        
        delta_xi_double   = norm(full( d_star),inf); 
        delta_dual_double = norm(full( (dual_k1 - dual_k)),inf);
        
        
          printf(obj.log,'debug','%-15d%-15e%-20e%-20e%-20e%-15.4f%-18e\n',...
                 i, new_cost , delta_xi_double, delta_dual_double, new_conVio , alpha_k , full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_k1,p0))))    );
        
        else
        
            delta_xi_double   = norm( full( d_star    ), inf ); 
            delta_dual_double = norm( full( dual_star ), inf );
            
          printf(obj.log,'debug','%-15d%-15e%-20e%-20e%-20e%-15.4f%-18e\n',...
                 i, new_cost, delta_xi_double, delta_dual_double, new_conVio, alpha_k , full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0))))   );
        
        end % --- end of display output ---
            
        %% check convergence criteria
        optimality_flag =     max ([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_k1,p0)))),new_conVio]) <= 1e-5;
        
        % check if solution stays below tolerance for a certain number of
        % iterations --> solved to acceptable level
        if max ([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_k1,p0)))),new_conVio]) <= 1e-4
            counter = counter + 1;
        else
            counter = 0;
        end % --- end of counter ---

        
        %% termination criteria and output preparation
        if i ~= obj.opts.max_iter && optimality_flag
                 
                  % store iteration info
                  info(i+1:end) = [];
                  obj.info.iter = info;
               
                   % adjust last solution similar to iteration i.e. overwrite optimization solution 
                   sol{1} = xi_k1;
                   sol{5} = dual_k1;
               
                   argout = sol;
        
                   printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
                   printf(obj.log,'debug','Solution status: Optimal solution found\n'); 
                   printf(obj.log,'debug',['Solution time: ' num2str(toc) ' s\n']); 
        
                   % terminate
                   return

        elseif i ~= obj.opts.max_iter && counter == 10
    
                           % store iteration info
                  info(i+1:end) = [];
                  obj.info.iter = info;
               
                   % adjust last solution similar to iteration i.e. overwrite optimization solution 
                   sol{1} = xi_k1;
                   sol{5} = dual_k1;
               
                   argout = sol;
        
                   printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
                   printf(obj.log,'debug','Solution status: Solved to acceptable level\n'); 
                   printf(obj.log,'debug',['Solution time: ' num2str(toc) ' s\n']); 
        
                   % terminate
                   return
        end % --- end of output preparation
       
        
        %% Damped BFGS 
        
        % update BFGS
        s = d_star; 
        y = casadi.DM( full(obj.nabla_xi_L(xi_k1,dual_k1,p0) - obj.nabla_xi_L(xi_k,dual_k1,p0)))';
           
        B = dampedBFGS(obj,B,s,y);
        
        
        %% new solution becomes new linearization point
        xi_k   = xi_k1;
        dual_k = dual_k1;

    end % --- end for-loop ---

    argout = sol;

end % --- end of function ---
