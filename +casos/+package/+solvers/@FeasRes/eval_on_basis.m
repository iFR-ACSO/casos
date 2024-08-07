function argout = eval_on_basis(obj,argin)

    % Evaluate sequential SOS problem.
    import casos.package.UnifiedReturnStatus
    
    args = argin;

      
    % initialize iteration info
    info = cell(1,obj.opts.max_iter);
    
    % Prepare display output
    printf(obj.log,'debug','%-15s%-15s%-20s%-20s%-20s%-15s%-10s\n', 'Iteration', 'cost', '||pr||_inf', '||du||_inf', 'squared_L2_vio','alpha','||dLdx||');
    printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');

    % first initial guess must be provided by user! 
    xi_k = args{1};

    % parameter from nonlinear problem
    p0   = args{2};
    
    % initialize Hessian/BFGS approximation
    B = eye(obj.sizeHessian(1));
 
    dual_k = [];
    
   conVio_xi0   = 0;
    
   % initialize the current cost with the last constraint violation of the
   % original problem
   curr_cost   = full(casadi.DM(p0(2)));
   curr_conVio = full(casadi.DM(p0(2)));
    
   % initialize filter
   obj.Filter =  obj.Filter.initializeFilter(conVio_xi0,curr_cost);

   counter = 0;

   tic
   % solve nonconvex SOS problem via sequence of convex SOS problem
    for i = 1:obj.opts.max_iter

        % initial guesss; actually not used but just for completness
        args{1} = xi_k;
    
        % set parameter to convex problem
        args{2}  = [p0; xi_k; B(:)];

        % adjust bounds
        args{3}  = argin{3}-xi_k;
        args{4}  = argin{4}-xi_k;

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

                   error('Feasibility restoration bubproblem infeasible! ')

            otherwise   % problem status unknown or ill-posed

                    d_star    = sol{1};
                    dual_star = sol{5};
        end


        %% backtracking line search filter 
        alpha_max            = 1;
        tau                  = 0.5;
        
        alpha_k = alpha_max;
		
        while true

            % compute alpha min 
            dirVecDeriv = full(casadi.DM(full(obj.nabla_xi_f(xi_k,p0)))*d_star);
            if dirVecDeriv < 0
                s_phi = 2.3;
                s_theta =  1.1;
                gamma_phi = 1;
                delta = 1;
                alpha_min = min([ 1e-5, ( gamma_phi*full(casadi.DM(full(obj.f(xi_k,p0)))) )/abs(dirVecDeriv ), (delta*full(casadi.DM(full(obj.f(xi_k,p0))))^s_theta )/(abs(dirVecDeriv)^s_phi) ]);

            else

                alpha_min = 1e-5; % same value as gamma_theta in filter

            end
                
            % leave linesearch 
            if alpha_k < alpha_min
                infeasible = 1;
                break
            else
                infeasible = 0;
            end


            % new-iterate
            xi_k1 = xi_k + alpha_k*d_star;

            % check cost, gradient of cost and constraint violation
            new_cost     = full(casadi.DM(full(obj.f(xi_k1,p0))));
            nabla_xi_f   = full(casadi.DM(full(obj.nabla_xi_f(xi_k1,p0))))';

            % measTime_Proj_in = tic;
            % xi_k1_sub = obj.xk1fun(xi_k1);
            % solPara_proj = obj.projConPara('p',[xi_k1_sub(1:4)]);
            % new_conVio   = sqrt(full(solPara_proj.f));
            % info{i+1}.filter_stats.measTime_proj_out = toc(measTime_Proj_in);

            new_conVio = 0 ;

            % check filter acceptance
            [obj.Filter,~,accept_new_iter] = obj.Filter.updateFilter(alpha_k,...
                                                                     full(d_star), ...
                                                                     curr_cost, ....
                                                                     curr_conVio,....
                                                                     new_cost,...
                                                                     new_conVio,....
                                                                     nabla_xi_f );

            
            if accept_new_iter
                break
            else
                % not acceptable to filter, adjust step length
                alpha_k = alpha_k*tau ;
            end
        
        
        end

        % compute dual variables
        if ~isempty(dual_k)
            dual_k1 = dual_k + alpha_k*(dual_star - dual_k);
        else
            dual_k1 = dual_star;
        end
                
        if infeasible
        
            

        end

        % set new iterate for next iteration
        curr_conVio = new_conVio;
        curr_cost   = new_cost;
        
        %% Prepare display output 
        if ~isempty(dual_k)
        
        delta_xi_double   = norm(full(casadi.DM(full( xi_k1 - xi_k))),inf); 
        delta_dual_double = norm(full(casadi.DM(full( dual_k1 - dual_k))),inf);
        
        
          printf(obj.log,'debug','%-15d%-15e%-20e%-20e%-20e%-15.4f%-18e\n',...
                 i, new_cost , delta_xi_double, delta_dual_double, new_conVio , alpha_k , full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k1,dual_k1,p0))))    );
        
        else
        
            delta_xi_double   = norm( full( d_star    ), inf ); 
            delta_dual_double = norm( full( dual_star ), inf );
            
          printf(obj.log,'debug','%-15d%-15e%-20e%-20e%-20e%-15.4f%-18e\n',...
                 i, new_cost, delta_xi_double, delta_dual_double, new_conVio, alpha_k , full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0))))   );
        
        end
            
        %% check convergence criteria
        % optimality_flag = new_cost/full(casadi.DM(p0(2))) < 1e-2;

        optimality_flag =   max([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0)))), new_cost]) <= 1e-3;
        
        % check if solution stays below tolerance for a certain number of
        % iterations --> solved to acceptable level
        if max([full(casadi.DM(full(obj.nabla_xi_L_norm(xi_k,dual_star,p0)))), new_cost]) <= 1e-2
            counter = counter + 1;
        else
            counter = 0;
        end

        
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

        elseif i == obj.opts.max_iter
                
                             % store iteration info
                  info(i+1:end) = [];
                  obj.info.iter = info;
               
                   % adjust last solution similar to iteration i.e. overwrite optimization solution 
                   sol{1} = xi_k1;
                   sol{5} = dual_k1;
               
                   argout = sol;
        
                   printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
                   printf(obj.log,'debug','Solution status: Maximum number of iterations reached\n'); 
                   printf(obj.log,'debug',['Solution time: ' num2str(toc) ' s\n']); 
        
                   % terminate
                   return
        end
       
        
        %% Damped BFGS 
        
        % update BFGS
        s = d_star; 
        y = casadi.DM( full(obj.nabla_xi_L(xi_k1,dual_star,p0) - obj.nabla_xi_L(xi_k,dual_star,p0)))';
           
        B = dampedBFGS(obj,B,s,y);
        
        
        %% new solution becomes new linearization point
        xi_k   = xi_k1;
        dual_k = dual_k1;

    end % --- end for-loop ---

    argout = sol;

end
