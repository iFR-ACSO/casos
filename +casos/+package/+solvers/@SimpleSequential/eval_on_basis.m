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
printf(obj.log,'debug','%-8s%-15s%-15s%-15s%-15s%-10s%-10s\n', 'iter', 'obj', ['||' char(916) 'x'  '||_inf'], ['||' char(916) char(955)  '||_inf'], [char(920) '(x_k)'],char(945),'||dLdx||');
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');

% input arguments
args = argin;

% initial guess (user)
x_k = args{1};

% parameter from nonlinear problem
p0   = args{2};

% initialize iteration info struct
info = cell(1,obj.opts.max_iter);

BFGS_scaling = 1;

% initialize BFGS-matrix
Bk     = eye(obj.init_para.size_B)*BFGS_scaling;

dual_k = zeros(obj.init_para.no_dual_var,1);

H = full(obj.hess_fun(x_k,p0,dual_k));

[L,D] = ldl(H);
D_PD = diag(max(eig(D),1e-8));

Bk0 = L*D_PD*L';

% initialize filter
args_conVio     =  args;
args_conVio{2}  =  [p0; x_k];
args_conVio{3}  = -inf(obj.init_para.conVio.no_con,1);
args_conVio{4}  =  inf(obj.init_para.conVio.no_con,1);

sol_convio = eval_on_basis(obj.solver_conVio, args_conVio);

theta_x_k = full(max(0,max(sol_convio{1})));
f_x_k     = full(obj.eval_cost(x_k,p0));

L_k = full(obj.L(x_k,p0,dual_k));
filter = [max(1,theta_x_k)*10, inf];

iter = 1;

% just for display output in first iteration
delta_xi_double   = norm(full( casadi.DM(x_k)),inf);
delta_dual_double = norm(full( (dual_k)),inf);
alpha_k           = 1;

feasibility_flag = 1;

eps_conVio = 1e-6; % about sqrt(eps)
eps_opt    = 1e-4;

measTime_seqSOS_in = tic;
while iter <= obj.opts.max_iter
    
    % display output 
    if ~mod(iter,10) && iter > 0
        printf(obj.log,'debug','%-8s%-15s%-15s%-15s%-15s%-10s%-10s\n', 'iter', 'obj', ['||' char(916) 'x'  '||_inf'], ['||' char(916) char(955)  '||_inf'], [char(920) '(x_k)'],char(945),'||dLdx||');
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
    end
    

     % scale = max(100,norm(full(dual_k),1)/length(dual_k))/100;

    printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
            iter-1, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , full(casadi.DM( full(obj.eval_gradLang(x_k,p0,dual_k)) ))   );
        
    %% check convergence (first-order optimality)
    if full(obj.eval_gradLang(x_k,p0,dual_k))  <= eps_opt*max(1,full(obj.eval_gradLang(x_k,p0,dual_k))) && theta_x_k <= eps_conVio
       
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Optimal solution found\n');
        solveTime = toc(measTime_seqSOS_in);
        printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);
        feasibility_flag = 1;

        break

    end
    

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

                
    % Invoke feasibility restoration if necessary
    if feas_res_flag >= 1 && ... % restoration requested
       theta_x_k > eps_conVio    % is the constraint violation larger then tolerance; otherwise restoration does not make sense
        
        % augment filter with the iterate at which we invoke restoration
        gamma_theta = 1e-3;
        gamma_phi   = 1e-3;

        % augment filter 
        filter = [filter;[theta_x_k*(1-gamma_theta), f_x_k - gamma_phi*theta_x_k]];
        
        % feasibility restoration phase
        [restored_sol,~,filter,feasibility_flag] = obj.feas_res_solver.eval_extended(filter,x_k,theta_x_k,p0,obj);

        if feasibility_flag >=1
           % continue from restored solution
           delta_xi_double   = norm(full( restored_sol.x_k1 - x_k) ,inf);
           delta_dual_double = norm(full( restored_sol.dual_k1 - dual_k ),inf);

           x_k       = restored_sol.x_k1;
           dual_k    = restored_sol.dual_k1;
           theta_x_k = restored_sol.theta_x_k1;
           f_x_k     = restored_sol.f_x_k1;
           alpha_k   = restored_sol.alpha_k;
            
      
        else
            % could not recover return iterate where restoration was
            % invoked
            sol    = sol_qp;
        
            sol{1} = x_k;
            sol{5} = dual_k;

            break
        end
        
    elseif feas_res_flag >= 1 &&  theta_x_k < eps_conVio
            % we are already better then treshhold i.e. restoration can't
            % make it better
            sol    = sol_qp;

            % if one, means Q-SDP is infeasible, thus sol.iter should be
            % zero hence return "old" iterate
            if feas_res_flag == 1
                sol{1} = x_k;
                sol{5} = dual_k;
            else
                 sol{1} = sol_iter.x_k1;
                sol{5} = sol_iter.dual_k1;
            end

            % stalled
            feasibility_flag = -2;

            break



    else

       % for display output
       delta_xi_double   = norm(full( sol_iter.x_k1 - x_k) ,inf);
       delta_dual_double = norm(full( sol_iter.dual_qp - dual_k ),inf);

       x_k       = sol_iter.x_k1;
       dual_k    = sol_iter.dual_qp;
       theta_x_k = sol_iter.theta_x_k1;
       f_x_k     = sol_iter.f_x_k1;
       alpha_k   = sol_iter.alpha_k;
    

    end

    % return last solution
    sol    = sol_qp;
    
    % exchange the primal/dual solution of QP with the accepted iterate
    sol{1} = sol_iter.x_k1;
    sol{5} = sol_iter.dual_k1;


    iter = iter + 1;

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

if feasibility_flag == 0
        
    % print last iterate
        printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
            iter, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , full(casadi.DM( full(obj.eval_gradLang(x_k,p0,dual_k)) ) )   );
        
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Problem infeasible.\n');
        solveTime = toc(measTime_seqSOS_in);
        printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);  

elseif feasibility_flag == -1


        % print last iterate
        printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
            iter, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , full(casadi.DM( full(obj.eval_gradLang(x_k,p0,dual_k)) ) )   );
        

        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Solver stalled in feasibility restoration.\n');
        solveTime = toc(measTime_seqSOS_in);
        printf(obj.log,'debug',['Solution time: ' num2str(solveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+solveTime) ' s\n']);

elseif feasibility_flag == -2


        % print last iterate
        printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
            iter, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , full(casadi.DM( full(obj.eval_gradLang(x_k,p0,dual_k)) ) )   );
        

        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Feasible solution. Restoration not invoked because already feasible.\n');
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
