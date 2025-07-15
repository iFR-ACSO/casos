function [sol_iter,sol,filter_glob,feasibility_flag,info_glob] = eval_extended(obj,filter_glob,x_R,theta_x_R0,p0,solver,info_glob,iter_glob)

% original parameter
p00 = p0;


% parameter for feasibility restoration: actual parameter and iterate where
% the restoration was/is invoked
p0 = [p00; x_R;0];

% input arguments
args = cell(10,1);

args{1} = x_R;
args{2} = x_R;
args{3} = -inf(obj.sparsity_pat_para.sparsity_xl.nnz,1);
args{4} =  inf(obj.sparsity_pat_para.sparsity_xl.nnz,1);
args{6} = -inf(obj.sparsity_pat_para.sparsity_gl.nnz,1);
args{7} =  inf(obj.sparsity_pat_para.sparsity_gl.nnz,1);
     
% initial guess for 
r0 = ones(solver.FeasRes_para.n_r,1);

% decision variabbles of feasibility restoration
x_k = [r0; x_R];

% initialize iteration info struct
info = cell(1,obj.opts.max_iter);


dual_k = zeros(obj.init_para.no_dual_var,1);

% initialize BFGS-matrix for feasibility restoration
% if strcmp(obj.opts.Hessian_init,'Identity')

    Bk =  eye(obj.init_para.size_B)*obj.opts.scale_BFGS0;

% elseif strcmp(obj.opts.Hessian_init,'Analytical')
% 
%     
% 
%     H = full(obj.hess_fun(x_k,p0,dual_k));
% 
%     Bk = casos.package.solvers.SequentialCommon.regularizeHessian(H);
% 
% end


% initialize filter for feasibility restoration
args_conVio     =  cell(10,1);
args_conVio{2}  =  [p0;  x_k];
args_conVio{3}  = -inf(obj.init_para.conVio.no_con,1);
args_conVio{4}  =  inf(obj.init_para.conVio.no_con,1);

sol_convio = eval_on_basis(obj.solver_conVio, args_conVio);

theta_x_k = full(max(0,max(sol_convio{1})));
f_x_k     = full(obj.eval_cost(x_k,p0));


filter = [max(1,theta_x_k)*10, inf];

lambda_min = 0.01;
lambda_max = 1;
lambda0 = lambda_min + (lambda_max-lambda_min)*1/(1+(max(theta_x_k,obj.opts.tolerance_con)-obj.opts.tolerance_con)/obj.opts.tolerance_con);

% parameter for feasibility restoration: actual parameter and iterate where
% the restoration was/is invoked
p0 = [p00; x_R;lambda0 ];



iter      = 1;
kappa_res = 0.9;

% just for display output in first iteration
delta_xi_double   = norm(full( casadi.DM(x_k)),inf);
delta_dual_double = norm(full( (dual_k)),inf);
alpha_k            = 1;

filter_Acceptance = 0;
feasibility_flag = 0;

measTime_seqSOS_in = tic;
while iter <= obj.opts.max_iter
    
    printf(obj.log,'debug','r%-8dr%-15er%-15er%-15er%-15er%-10fr%-10e\n',...
            iter-1, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , full(casadi.DM( full(obj.eval_gradLang(x_k,p0,dual_k)) ) )   );

        
    %% check convergence (first-order optimality) of filter
    if full(obj.eval_gradLang(x_k,p0,dual_k))  <= obj.opts.tolerance_opt*max(1,full(obj.eval_gradLang(x_k,p0,dual_k))) ...
        && theta_x_k <= obj.opts.tolerance_con && feasibility_flag == 0
       
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Optimal solution found in feasibility restoration. Problem seems infeasible.\n');
        solveTime = toc(measTime_seqSOS_in);
        printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);
        
        % if the feasibility restoration found an optimal solution and not
        % acceptable to the filter of the original problem
        feasibility_flag = 0;

        break

    end
    
    feasibility_flag = 0;

    % compute solution of current iterate: solve Q-SDP, perform linesearch and update BFGS
    [sol_iter,sol_qp,feas_res_flag,info,obj,filter,Bk] = do_single_iteration(obj, ...
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

   %% check feasibility; 
   % feasibility restoration does not have a feasibility restoration
   if feas_res_flag == 1
       % QP  of feasibility restoration is also infeasible --> Problem
       % seems infeasible
       feasibility_flag = 0;
       break
   elseif feas_res_flag == 2
       % alpha below alpha min i.e. solver stalled
       feasibility_flag = -1;
       break
   else 

       %% check acceptance to filter of acutal problem

      % extract solution from feasibility restoration
      x_k       = sol_iter.x_k1;
      dual_k    = sol_iter.dual_k1;
      alpha_k   = sol_iter.alpha_k;

      % first n_r decision variables are the individual constraint violations      
      theta_x_k_or = full(max(0,max(x_k(1:solver.FeasRes_para.n_r))));
      
      % store them
      info{iter}.constraint_violation = x_k(1:solver.FeasRes_para.n_r);

      % cost (original problem) at trial point 
      f_x_k1_or   = full(solver.eval_cost(x_k(solver.FeasRes_para.n_r+1:end),p00));
        
     % check acceptance to filter of the original problem
     theta_l = filter_glob(:,1);
     f_l     = filter_glob(:,2);
            
     % new point lies in forbidden region if both are larger than filter entries
     dominance_bool      = [];
     dominance_bool(:,1) = theta_x_k_or >= theta_l; 
     dominance_bool(:,2) = f_x_k1_or     >= f_l;
            
     % check pairs; if one means both lie in forbidden region
     dominance_bool = all(dominance_bool, 2);
            
     if any(dominance_bool) % means not acceptable to filter
         filter_Acceptance = 0;
     else
         filter_Acceptance = 1;
     end
     
    
    % check acceptance to filter of original problem plus threshold for constraint violation
    if filter_Acceptance && ...              % acceptable to filter          
        theta_x_k_or <= theta_x_R0*kappa_res % threshold

        feasibility_flag = 1;
        break

    end

   end

    % lambda0 = lambda_min+(lambda_max-lambda_min)*1/(1+ (max(theta_x_k_or, obj.opts.tolerance_con) - obj.opts.tolerance_con)/obj.opts.tolerance_con);

    % parameter for feasibility restoration: actual parameter and iterate where
    % the restoration was/is invoked
    p0 = [p00; x_R;lambda0];


       % for display output
       delta_xi_double   = norm(full( sol_iter.x_k1 - x_k) ,inf);
       delta_dual_double = norm(full( sol_iter.dual_qp - dual_k ),inf);

       x_k       = sol_iter.x_k1;
       dual_k    = sol_iter.dual_k1;
       theta_x_k = sol_iter.theta_x_k1;
       f_x_k     = sol_iter.f_x_k1;
       alpha_k   = sol_iter.alpha_k;
    

   iter = iter + 1;

end % end of while


if iter >= obj.opts.max_iter
    % if max. number of iterations reached, the solver stalled
    feasibility_flag = -1;
end

if feasibility_flag <= 0
            
   % return last solution of q-SDP
   sol    = sol_qp;
        
   sol{1} = x_k(solver.FeasRes_para.n_r+1:end);
   sol{5} = dual_k(solver.FeasRes_para.n_r+1:end);

else
    
    % return solution of successful restoration
    sol    = sol_qp;
    
    sol{1} = casadi.DM(x_k(solver.FeasRes_para.n_r+1:end));
    sol{5} = dual_k(solver.FeasRes_para.length_dualOut+1:end);

    % output to original problem
    sol_iter.x_k1       = sol_iter.x_k1(solver.FeasRes_para.n_r+1:end);
    sol_iter.dual_k1    = sol_iter.dual_k1(solver.FeasRes_para.length_dualOut+1:end);
    sol_iter.f_x_k1     = f_x_k1_or ;
    sol_iter.theta_x_k1 = theta_x_k_or;
   
end

info_glob{iter_glob} = info{iter};

end % end of eval_extended

