function argout = call(obj,argin)
    % Call sequential sum-of-squares
    
    import casos.package.UnifiedReturnStatus
    
    % get arguments
    args = argin;
    
    % initialize iteration info
    info = cell(1,obj.opts.max_iter);
    
    % first initial guess must be provided by user!
    xi_k = args{1};

    % check if we have parameter from the original nonlinear problem
    if ~is_zero(args{2})
        p0   = args{2};
    else
        p0 = [];
    end
    
    if obj.opts.verbose 
        % disp('Start sequential SOS...')
        printf(obj.log,'debug','%-15s%-15s%-20s%-20s%-20s%-20s\n', 'Iteration', 'cost', 'squared_L2_pr', 'squared_L2_du', 'squared_L2_vio','alpha');
        printf(obj.log,'debug','-----------------------------------------------------------------------------------------------\n');
    end

    % solve nonconvex SOS problem via sequence of convex SOS problem
    for i = 1:obj.opts.max_iter
        
        % initial guesss; actually not used but just for completness
        args{1} = xi_k;
    
        % set parameter to convex problem
        % args{2}  = [p0; xi_k];
        args{2} = obj.vertcatFun(p0, xi_k);
    
        % evaluate convex SOS problem
        sol = call(obj.sossolver, args);
  
        % store iteration info
        info{i} = obj.sossolver.stats;
    
    
        %% check solution status from conic optimization
        switch (info{i}.UNIFIED_RETURN_STATUS)
            case UnifiedReturnStatus.SOLVER_RET_SUCCESS
                
                % set initial guess and finally parameter; 
                % sos constraint reads  delta x = f(x0) + df/dx(x0)*(x-x0);
                % where delta x = xi_k+1 - xi_k hence xi_k+1 = delta_x + xi_k 

                % xi_plus   = sol{1} + xi_k;
                xi_plus = obj.plusFun(sol{1},xi_k);

                dual_plus = sol{5};
                
                % constraint_vio = sqrt(full(casadi.DM( pnorm2( obj.constraintFun(xi_plus,p0) - sol{3} ) )));
                 
                % delta_constraint = obj.constraintFun(xi_plus,p0) - sol{3} ;
                % constraint_vio  = full(casadi.DM( dot( delta_constraint, delta_constraint ) ));
                % constraint_vio  = full(casadi.DM(obj.norm2FunVio(xi_plus,p0,sol{3})));
                
                sos_g = sol{3};
                
                % solve projection problem to get constraint violation
                % try 
                sol_proj = obj.projCon('p',sos_g(obj.idxNonlinCon));
          
                constraint_vio  = double(sol_proj.f);
                % catch
                    % constraint_vio  = full(casadi.DM(obj.norm2FunVio(xi_plus,p0,sol{3})));
                % end
                cost            = double(sol{2});

            case UnifiedReturnStatus.SOLVER_RET_UNKNOWN
                 % problem seems feasible, but but solution is not

                xi_plus   = sol{1} + xi_k;
    
                dual_plus = sol{5};
                
                % constraint_vio = sqrt(full(casadi.DM( pnorm2( obj.constraintFun(xi_plus,p0) - sol{3} ) )));
                 
                % delta_constraint = obj.constraintFun(xi_plus,p0) - sol{3} ;
                constraint_vio  = full(casadi.DM( dot( delta_constraint, delta_constraint ) ));
                cost            = double(sol{2});
                
            case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE

                 % if underlying conic problem is infeasible go to restoration phase   
                 if i > 1
    
                   % store iteration info
                   info(i+1:end) = [];
                   obj.info.iter = info;
                    
                   % take solution from previous iteration
                   argout = sol_old;

                    printf(obj.log,'error','Convex optimization is infeasible.');
                    assert(~obj.opts.error_on_fail,'Convex optimization run into numerical errors.')
    
                  return
    
                else
    
                    % error: failed
                    obj.status = UnifiedReturnStatus.SOLVER_RET_NAN;
                    assert(~obj.opts.error_on_fail,'Problem is primal and/or dual infeasible or unbounded.')   
    
                end

            otherwise   % problem status unknown or ill-posed
                
    
                if i > 1
    
                   % store iteration info
                   info(i+1:end) = [];
                   obj.info.iter = info;
                    
                   % take solution from previous iteration
                   argout = sol_old;

                 printf(obj.log,'error','Convex optimization run into numerical errors.');
                 assert(~obj.opts.error_on_fail,'Convex optimization run into numerical errors.')
    
                  return
    
                else
    
                    % error: failed
                    obj.status = UnifiedReturnStatus.SOLVER_RET_NAN;
                    assert(~obj.opts.error_on_fail,'Convex optimization run into numerical errors.')
    
                end
      
            
        end
    

        %% Filter
        if i == 1
           % initialize
           Filter = [cost , constraint_vio];
                
        else
           
       
           List_cost   = Filter(:,1);
           List_conVia = Filter(:,2);
                
                % check if any point in the list dominates the new iterate
                if any(List_cost < cost) && any(List_conVia < constraint_vio)
                    FilterAccept = 0;
                else
                    FilterAccept = 1;
                    Filter(i,1) = cost ;
                    Filter(i,2) = constraint_vio;
                end
        end
    
        
        if i == 1
             % set current solution as previous solution for next iteration
             xi_k    = xi_plus;
             dual_k  = dual_plus;
             sol_old = sol;
             continue
        end


               % perform linesearch if curren iterate is not accepted to
               % filter
               if ~FilterAccept
                    % select an algorithm to perform the linesearch
                    switch(obj.opts.line_search) 
                        case 'fminbnd'
                                dopt = line_search_fminbnd(obj, ...
                                                           p0,...
                                                           xi_k, ...
                                                           xi_plus, ...
                                                           dual_plus );
                        case 'bisection'
                                dopt = bisection_minimization(obj, ...
                                                              p0,...
                                                              xi_k, ...
                                                              xi_plus, ...
                                                              dual_plus );
                        case 'polySolver'
                                sol = obj.lineSearch('p',[xi_k; xi_plus; dual_plus],...
                                                     'lbx',0.1,...
                                                     'ubx',1);
    
                                dopt = double(sol.x);
    
                    end
               else

                     dopt = 1;
                end

                % update primal and dual solution
                [xi_k1, dual_k1] = obj.updateLineSearch(xi_k, xi_plus, dual_k, dual_plus, dopt );
                % xi_k1    = dopt*xi_plus    + (1-dopt)*xi_k;
                % duals_k1 = dopt*dual_plus  + (1-dopt)*dual_k;
                 
            
            % check convergence or if maximum number of iterations are reached    
            [delta_xi, delta_dual] = obj.deltaOptVar(xi_k,xi_k1,dual_k,dual_k1);
            % delta_xi    = xi_k1    - xi_k;
            % delta_dual  = dual_k1 - dual_k;

            % delta_xi_double   = full( casadi.DM( (dot(delta_xi,delta_xi)) ) );
            % delta_dual_double = full( casadi.DM( (dot(delta_dual,delta_dual))) );
            
            [delta_xi,delta_dual] = obj.norm2FunOptVar(delta_xi,delta_xi,delta_dual,delta_dual);
            delta_xi_double   = full( casadi.DM(delta_xi));
            delta_dual_double = full( casadi.DM(delta_dual) );

             if obj.opts.verbose 
                printf(obj.log,'debug','%-15d%-15f%-20e%-20e%-20e%-20.4f\n',...
                        i, cost , delta_xi_double, delta_dual_double, constraint_vio, dopt  );
            end
    

            if i == obj.opts.max_iter || ...
               delta_xi_double      < obj.opts.tolerance_abs && ...
               delta_dual_double    < obj.opts.tolerance_rel*full( casadi.DM(  (dot(dual_k1,dual_k1)) ) )  
               % full( casadi.DM( pnorm2()) )       < obj.opts.tolerance_abs && ...
               % full( casadi.DM( pnorm2(duals_k1 - dual_k) ) ) < obj.opts.tolerance_rel*full( casadi.DM( pnorm2( duals_k1 ) ) ) 
              
               % store iteration info
               info(i+1:end) = [];
               obj.info.iter = info;
        
                % adjust last solution similar to iteration i.e. overwrite optimization solution 
                sol{1} = xi_k1;
                sol{5} = dual_k1;
        
                argout = sol;
                printf(obj.log,'debug','---------------------------------------------------------------------------------\n');
                printf(obj.log,'debug','Solution status: Optimal solution found\n'); 

                % terminate
                return
        
            end % end-if convergence check

            % set current solution as previous solution for next iteration
            xi_k      = xi_k1;
            dual_k    = dual_k1;
            sol_old   = sol;
    
    
    end % end for-loop sequential sos
end % end of function
