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
    
    dual_k = args{end};
    
    % check if an initial guess for duals is available
    if is_zero(dual_k)
        dual_k = [];
    end

     % start sequential SOS
    printf(obj.log,'debug','%-15s%-15s%-20s%-20s%-20s%-15s%-10s\n', 'Iteration', 'cost', 'squared_L2_pr', 'squared_L2_du', 'squared_L2_vio','alpha','||dLdx||');
    printf(obj.log,'debug','----------------------------------------------------------------------------------------------------------------------\n');

    SDP_inf_flag = 0;

    % initilize Hessian if SQP method is used
    if strcmp(obj.opts.Sequential_Algorithm,'SLP')
        B = [];
    else
         B = eye(obj.sizeHess(1));
    end
    

    % initialize filter using the initial guess
    solPara_proj = obj.projConPara('p',[xi_k;p0]);

    con_vio_0 = (double(solPara_proj.f));
    cost_0    = inf;

    Filter = [cost_0 , con_vio_0];
        
    oldcost    = cost_0;
    old_conVio = con_vio_0;
    
    tic
    % solve nonconvex SOS problem via sequence of convex SOS problem
    for i = 1:obj.opts.max_iter

        stats.solver_stats      = [];
        stats.iter              = 0;
        stats.dual              = [];
        stats.primal            = [];
        stats.cost              = [];
        stats.conVio            = 0;
        stats.nabla_Langrange   = [];
        
        % initial guesss; actually not used but just for completness
        args{1} = xi_k;
    
        % set parameter to convex problem
        
        args{2}  = [p0; xi_k; B(:)];


        % adjust bounds
        args{3}  = argin{3}-xi_k;
        args{4}  = argin{4}-xi_k;
        % args{5}  = argin{5}-xi_k(end);

        % evaluate convex SOS problem
        sol = call(obj.sossolver, args);
  
        % store iteration info
        stats.solver_stats = obj.sossolver.stats;
        % info{i} = obj.sossolver.stats;
        
        
        %% check solution status from conic optimization
        switch (obj.sossolver.stats.UNIFIED_RETURN_STATUS)
            case UnifiedReturnStatus.SOLVER_RET_SUCCESS
                
                % set initial guess and finally parameter; 
                % sos constraint reads  delta x = f(x0) + df/dx(x0)*(x-x0);
                % where delta x = xi_k+1 - xi_k hence xi_k+1 = delta_x + xi_k 

                % xi_plus   = sol{1} + xi_k;
                xi_plus = obj.plusFun(sol{1},xi_k);
                
                dual_plus = sol{5};

            case UnifiedReturnStatus.SOLVER_RET_UNKNOWN
                 % problem seems feasible, but but solution is not

                xi_plus   = sol{1} + xi_k;
    
                dual_plus = sol{5};
                
                
            case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE
                
                 SDP_inf_flag = 1;

            otherwise   % problem status unknown or ill-posed
                
                xi_plus = obj.plusFun(sol{1},xi_k);
                
                dual_plus = sol{5};
            
        end
    
    if ~SDP_inf_flag
       switch(obj.opts.globalization) 
           
           % naive implementation of simple linesearch with a filter
           case 'filter_linesearch_simple'
            
            [xi_k1, dual_k1, Filter, Linesearch_stall,alpha,predictedcost,predictedConVio] = ...
                lineSearchFilterSimple(obj,p0,Filter,xi_k,xi_plus,dual_k,dual_plus); 

           case 'filter_linesearch'
                

               
               [xi_k1, dual_k1, Filter, Linesearch_stall,alpha,predictedcost,predictedConVio]   = lineSearchFilter(obj,p0,Filter,xi_k,xi_plus,dual_k,dual_plus,old_conVio,B,args) ;  
               old_conVio = predictedConVio;
           
                    
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

        if ~isempty(dual_k)

            [delta_xi,delta_dual] = obj.norm2FunOptVar(xi_k,xi_k1,dual_k,dual_k1);

            delta_xi_double   = (full( casadi.DM(delta_xi))); 
            delta_dual_double = (full( casadi.DM(delta_dual) ));


          printf(obj.log,'debug','%-15d%-15f%-20e%-20e%-20e%-15.4f%-18e\n',...
                 i, predictedcost , delta_xi_double, delta_dual_double, predictedConVio , alpha , norm(double(obj.dLdx(xi_k1,dual_k1,p0)),inf)  );

        else

           [delta_xi,~] = obj.norm2FunOptVar(xi_k,xi_k1,dual_k1,dual_k1);

            delta_xi_double   = (full( casadi.DM(delta_xi))); 
            delta_dual_double = (full( casadi.DM(pnorm2(dual_k1)) ));


            printf(obj.log,'debug','%-15d%-15f%-20e%-20e%-20e%-15.4f%-18e\n',...
                i, predictedcost , delta_xi_double, delta_dual_double, predictedConVio , alpha , norm(double(obj.dLdx(xi_k1,dual_k1,p0)),inf)  );

        end
    

        optimality_flag = ~Linesearch_stall &&    norm(double(obj.dLdx(xi_k1,dual_k1,p0)),inf)   <= obj.opts.optimality_tol && ...
                                                  predictedConVio                                <= obj.opts.conVio_tol && ...
                                                  delta_dual_double                              <= obj.opts.dual_tol;
                                   
                                                  

            converged_flag = 0;                                                 


    end

        stats.iter              = i;
        stats.dual              = delta_dual_double;
        stats.primal            = delta_xi_double;
        stats.cost              = predictedcost;
        stats.nabla_Langrange   = norm(double(obj.dLdx(xi_k1,dual_k1,p0)),inf);
        stats.conVio            = predictedConVio;
        
       info{i} = stats;
        
       if SDP_inf_flag
                    
                   % store iteration info
                   info(i+1:end) = [];
                   obj.info.iter = info;
                   

                   % if UnifiedReturnStatus.SOLVER_RET_INFEASIBLE
                       printf(obj.log,'debug','Problem is primal and/or dual infeasible or unbounded.\n'); 
                       printf(obj.log,'debug',['Solution time: ' num2str(toc) ' s\n']);
                    
                  argout = sol;     
                  return

       elseif i ~= obj.opts.max_iter && optimality_flag
          
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

        elseif i ~= obj.opts.max_iter &&   Linesearch_stall

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
   
    
 end % end for-loop sequential sos
end % end of function
