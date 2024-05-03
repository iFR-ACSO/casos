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
    

     % start sequential SOS
    printf(obj.log,'debug','%-15s%-15s%-20s%-20s%-20s%-15s%-10s\n', 'Iteration', 'cost', 'squared_L2_pr', 'squared_L2_du', 'squared_L2_vio','alpha','||dLdx||');
    printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');

    SDP_inf_flag = 0;
    counter = 0;
    
    if strcmp(obj.opts.Sequential_Algorithm,'SLP')
        B = [];
    else
         B = eye(obj.sizeHess(1));
    end
   
    tic
    % solve nonconvex SOS problem via sequence of convex SOS problem
    for i = 1:obj.opts.max_iter
        
        % initial guesss; actually not used but just for completness
        args{1} = xi_k;
    
        % set parameter to convex problem
        
        args{2}  = [p0; xi_k; B(:)];
        % args{2} = obj.vertcatFun(p0, xi_k);

        % adjust bounds
        args{3}  = argin{3}-xi_k;
        args{4}  = argin{4}-xi_k;
        % args{5}  = argin{5}-xi_k(end);

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

                % norm(full(casadi.DM((poly2basis(sol{2} + dot(sol{5},sol{3}))))))
                % casadi.DM(obj.langrangeLinear(sol{1},xi_k,p0,dual_plus,B))
     

            case UnifiedReturnStatus.SOLVER_RET_UNKNOWN
                 % problem seems feasible, but but solution is not

                xi_plus   = sol{1} + xi_k;
    
                dual_plus = sol{5};
                
                sol_proj = obj.projCon('p',sos_g(obj.idxNonlinCon));
                
                constraint_vio  = double(sol_proj.f);
                cost            = double(sol{2});
                
            case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE
                


                 SDP_inf_flag = 1;


            otherwise   % problem status unknown or ill-posed
                
    
               SDP_inf_flag = 1;
            
        end
    

        %% Linesearch-Filter method 
        if i == 1
            
                sos_g     = sol{3};
                
                % solve projection problem to get constraint violation
                sol_proj = calcProj(obj,sos_g,obj.idxNonlinCon);
          
                predictedConVio  = double(sol_proj.f);
  
                predictedcost    = double(obj.cost_fun(xi_plus));

           % initialize filter in first iteration
            Filter = [predictedcost , predictedConVio];

             % set current solution as previous solution for next iteration
             xi_k    = xi_plus;
             dual_k  = dual_plus;
             oldcost   = predictedcost;
             sol_old = sol;
             
             continue
           
        end

       switch(obj.opts.globalization)
           
           % naive implementation of simple linesearch with a filter
           case 'filter_linesearch_simple'
            
            [xi_k1, dual_k1, Filter, Linesearch_stall,alpha,predictedcost,predictedConVio] = ...
                lineSearchFilterSimple(obj,p0,Filter,xi_k,xi_plus,dual_k,dual_plus); 

           case 'fminbnd'
              Linesearch_stall = 0;
              alpha =  line_search_fminbnd(obj, p0, xi_k, xi_plus, dual_plus );

            xi_k1   = xi_k   + alpha*(xi_plus - xi_k);
            dual_k1 = dual_k + alpha*(dual_plus - dual_k);

            solPara_proj = obj.projConPara('p',[xi_k1;p0]);

            predictedConVio = double(solPara_proj.f);

            predictedcost   = double(obj.cost_fun(xi_k1));
        
       end


        % Prepare display output    
        [delta_xi,delta_dual] = obj.norm2FunOptVar(xi_k,xi_k1,dual_k,dual_k1);

        delta_xi_double   = (full( casadi.DM(delta_xi))); 
        delta_dual_double = (full( casadi.DM(delta_dual) ));


      printf(obj.log,'debug','%-15d%-15f%-20e%-20e%-20e%-15.4f%-18e\n',...
             i, predictedcost , delta_xi_double, delta_dual_double, norm(predictedConVio,2) , alpha , norm(double(obj.dLdx(xi_k1,dual_k1,p0)),2)  );
    
        
        % compute convergence flags
        optimality_flag = ~Linesearch_stall && norm(predictedConVio)                      <= 1e-5 && ...
                                               norm(double(obj.dLdx(xi_k1,dual_k1,p0)))   <= 1e-3 && ...
                                               norm(predictedcost - oldcost)              <= 1e-6 && ...
                                               delta_xi_double                            <= 1e-6 && ...
                                               delta_dual_double                          <= 1e-6;


        if optimality_flag
            counter = 0;
        end

        converged_flag  =  ~Linesearch_stall  && norm(double(obj.dLdx(xi_k1,dual_k1,p0)))   > 1e-3 && ...
                                                             norm(predictedConVio)         <= 1e-5 && ...
                                                             norm(predictedcost-oldcost)   <= 1e-6 && ...
                                                             delta_xi_double               <= 1e-6 && ...
                                                             delta_dual_double             <= 1e-6;

        if converged_flag
            counter = counter +1;
        end


        stall_flag  =  Linesearch_stall  || norm(double(obj.dLdx(xi_k1,dual_k1,p0)))  > 1e-3 && ...
                                            norm(predictedConVio)                     > 1e-5 && ...
                                            norm(predictedcost-oldcost)               <= 1e-6 && ...
                                            delta_xi_double                           <= 1e-6 && ...
                                            delta_dual_double                         <= 1e-6;


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

        elseif i ~= obj.opts.max_iter && converged_flag 

                           % store iteration info
           info(i+1:end) = [];
           obj.info.iter = info;
        
            % adjust last solution similar to iteration i.e. overwrite optimization solution 
            sol{1} = xi_k1;
            sol{5} = dual_k1;
        
            argout = sol;
            printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
            printf(obj.log,'debug','Solution status: Solved to acceptable level!\n'); 
            printf(obj.log,'debug',['Solution time: ' num2str(toc) ' s\n']); 
            
            % terminate
            return

        elseif i ~= obj.opts.max_iter &&   stall_flag

                                       % store iteration info
           info(i+1:end) = [];
           obj.info.iter = info;
        
            % adjust last solution similar to iteration i.e. overwrite optimization solution 
            sol{1} = xi_k1;
            sol{5} = dual_k1;
        
            argout = sol;
            printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');
            
            if Linesearch_stall
                stall_string = 'Filter linesearch could not find a feasible solution ';
            elseif stall_flag && ~Linesearch_stall
                stall_string = 'could not find a feasible solution ';
            else
                stall_string = 'Filter linesearch could not find a feasible solution and no progess ';
            end


            printf(obj.log,'debug',['Solution status: Stalled because ' stall_string '\n']); 
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
            printf(obj.log,'debug','Solution status: Reached maximum number of iterations.\n'); 
            printf(obj.log,'debug',['Solution time: ' num2str(toc) ' s\n']); 

            % terminate
            return

        elseif SDP_inf_flag
                    
                   % store iteration info
                   info(i+1:end) = [];
                   obj.info.iter = info;
                    
                   % take solution from previous iteration
                   argout = sol_old;

                   % if UnifiedReturnStatus.SOLVER_RET_INFEASIBLE
                       printf(obj.log,'debug','Problem is primal and/or dual infeasible or unbounded.\n'); 
                       printf(obj.log,'debug',['Solution time: ' num2str(toc) ' s\n']); 
          
            return
   
        end % end-if convergence check
        
      
        % Damped BFGS; see Nocedal p. 536/537 and/or Biegler  page 44 and pages 102/103
        if strcmp(obj.opts.Sequential_Algorithm,'SQP')
            s = casadi.DM(casadi.SX(obj.delta_search(xi_k1,xi_k))); 
            y = casadi.DM(casadi.SX(obj.dLdx(xi_k1,dual_k1,p0) - obj.dLdx(xi_k,dual_k1,p0)));
    
            B = dampedBFGS(obj,B,s,y);
        end
       
        % set current solution as previous solution for next iteration
        xi_k      = xi_k1;
        dual_k    = dual_k1;
        oldcost   = predictedcost;
        sol_old   = sol;

    
 end % end for-loop sequential sos
end % end of function
