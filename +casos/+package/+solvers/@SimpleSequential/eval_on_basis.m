function argout = eval_on_basis(obj,argin)
% Solve sequential SOS problem.

import casos.package.UnifiedReturnStatus


%% print output with current problem and setting
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug',['\t\t Ca' char(931) 'oS-Nonlinear sum-of-squares optimization suite \n']);
printf(obj.log,'debug','\t\t Sequential Quadratic Sum-of-Squares Solver v1.0, December''24 \n');
printf(obj.log,'debug','\t\t GNU GENERAL PUBLIC LICENSE \n');
printf(obj.log,'debug','\t\t Institute of Flight Mechanics and Control, Univ. of Stuttgart \n');
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug','Problem:\n Decision variables = %d (coefficents) \n SOS constraints \t= %d \n Linear constraints = %d\n',obj.display_para.no_decVar   , obj.display_para.no_sosCon, 0 );
printf(obj.log,'debug','Settings:\n Solver:\t%s\n max_iter:\t%d \n Con. Violation Check: \t%s\n', obj.display_para.solver,obj.opts.max_iter,'signed-distance');
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug','%-8s%-15s%-15s%-15s%-15s%-10s%-10s\n', 'iter', 'obj', '||pr||_inf', '||du||_inf', '||conVio||_2','alpha','||dLdx||');
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');


% input arguments
args = argin;

% initial guess (user)
x_k = args{1};

% parameter from nonlinear problem
p0   = args{2};

% initialize iteration info struct
info = cell(1,obj.opts.max_iter);

% initialize BFGS-matrix
Bk     = eye(obj.init_para.size_B);
dual_k = ones(obj.init_para.no_dual_var,1);

% initialize filter
args_conVio     =  args;
args_conVio{2}  =  [p0; x_k];
args_conVio{3}  = -inf(obj.init_para.conVio.no_con,1);
args_conVio{4}  =  inf(obj.init_para.conVio.no_con,1);

sol_convio = eval_on_basis(obj.solver_conVio, args_conVio);

theta_x_k = full(max(0,max(sol_convio{1})));
f_x_k     = full(obj.eval_cost(x_k,p0));

filter = [theta_x_k, f_x_k];

iter = 1;

% just for display output in first iteration
delta_xi_double   = norm(full( casadi.DM(x_k)),inf);
delta_dual_double = norm(full( (dual_k)),inf);
alpha_k            = 1;


measTime_seqSOS_in = tic;
while iter <= obj.opts.max_iter
    
    % display output
    if ~mod(iter,10) && iter > 0
        printf(obj.log,'debug','%-8s%-15s%-15s%-15s%-15s%-10s%-10s\n', 'iter', 'obj', '||pr||_inf', '||du||_inf', '||conVio||_2','alpha','||dLdx||');
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');

    end
    
    printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
            iter, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , full(casadi.DM( full(obj.eval_gradLang(x_k,p0,dual_k)) ) )   );
        
        
    %% check convergence (first-order optimality)
    if full(obj.eval_gradLang(x_k,p0,dual_k)) <= 1e-5 && ...
       theta_x_k <= 1e-7
        
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Optimal solution found\n');
        solveTime = toc(measTime_seqSOS_in);
        printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);
       
        break

    end

    % compute solution of current iterate: solve Q-SDP and perform
    % linesearch
    [sol_iter,sol_qp,feas_res_flag,info,obj,filter] = do_single_iteration(obj, ...
                                                                          iter,...
                                                                          x_k,...
                                                                          dual_k,...
                                                                          theta_x_k,...
                                                                          f_x_k,...
                                                                          Bk,...
                                                                          p0,...
                                                                          args, ...
                                                                          filter, ...
                                                                          info);
                
    % Invoke feasibility restoration if necessary
    if feas_res_flag 
        % invoke feasibility restoration
        break
    else

       Bk = damped_BFGS(obj,Bk,x_k,p0,sol_iter);
    
       delta_xi_double   = norm(full( sol_iter.x_k1 - x_k) ,inf);
       delta_dual_double = norm(full( sol_iter.dual_k1 - dual_k ),inf);

       x_k       = sol_iter.x_k1;
       dual_k    = sol_iter.dual_k1;
       theta_x_k = sol_iter.theta_x_k1;
       f_x_k     = sol_iter.f_x_k1;
       alpha_k   = sol_iter.alpha_k;
    
       iter = iter + 1;
    end

    % return last solution
    sol =sol_qp;

    sol{1} = sol_iter.x_k1;
    sol{5} = sol_iter.dual_k1;

end % end of while loop

% display output for user
if iter >= obj.opts.max_iter
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Maximum number of iterations reached\n');
        solveTime = toc(measTime_seqSOS_in);
        printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);  
end


% assign output
argout = sol;

% store iteration info
obj.info.iter = info(1:iter-1);

end % end of function
