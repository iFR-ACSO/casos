function [sol_iter,sol_qp,feas_res_flag,info,obj,filter,Bk] = do_single_iteration(obj, ...
    iter,...
    x_k,...
    dual_k,...
    theta_xk,...
    f_xk,...
    Bk,...
    p0,...
    args, ...
    filter, ...
    info)

%% evaluate convex SOS problem and check feasibility
[x_star,dual_star,sol_iter,sol_qp,feas_res_flag,info] = solve_Q_SDP(obj,iter,x_k,p0,Bk,args,info);

% initilizing timing of backstepping components in case we have to switch
% to feasibility restoration
info{iter}.timeStats.HessApproxTime         = 0;
info{iter}.timeStats.totalBackstepTime      = 0;
info{iter}.timeStats.FilterAcceptTime       = 0;
info{iter}.timeStats.SuffDecreaseTime       = 0;
info{iter}.timeStats.SocTime                = 0;


if feas_res_flag
    return % leave single iteration and invoke feasibility restoration
end

%% backtracking filter-linesearch
alpha_max   = obj.opts.filter_struct.alpha_max;

gamma_theta = obj.opts.filter_struct.gamma_theta ;
gamma_phi   = obj.opts.filter_struct.gamma_phi ;

LangrangeFilter = obj.opts.filter_struct.LangrangeFilter;

alpha = alpha_max;

f_type = 0;
amijo  = 0;


measBackStepTime = tic;
SocTimeMeas      = [];
SuffDecreaseTimeMeas = 0;

% just for first iteration
alpha_min = 0;

% invoke feasibility restoration if step-length is below minimum
feas_res_flag = 2;
% Backtracking linesearch
while alpha > alpha_min 

    % compute search direction (primal and dual)
    dk   = (x_star    - x_k);
    dkl  = (dual_star-dual_k);


    measFilterAcceptTime = tic;
    % check filter acceptance
    [x_k1,theta_x_k1,f_x_k1 ,filter_Acceptance,sol_convio] = check_filter_acceptance(obj,filter,alpha,x_k,dual_k,dk,dkl,p0,args);
    FilterAcceptTimeMeas = toc(measFilterAcceptTime);

    dual_k1 = dual_k + alpha*(dual_star-dual_k);

    % heuristic: there might be blocking entries in the filter:
    % check first-order optimality conditions to potentially aboard
    % linesearch early as possible if we are "optimal"
    if full(obj.eval_gradLag(x_k1,p0,dual_k1))  <= obj.opts.tolerance_opt*max(1,full(obj.eval_gradLag(x_k1,p0,dual_k1))) && theta_x_k1 <= obj.opts.tolerance_con
        feas_res_flag = 0;
        break
    end

    if LangrangeFilter
        L_k1= full(obj.L(x_k1,p0,dual_k1));
        L_k = full(obj.L(x_k,p0,dual_k));
    else
        L_k1 = [];
        L_k  = [];
    end

    if filter_Acceptance

        measSuffDecreaseTime = tic;
        % check sufficient decrease
        [suffDecrease_flag,f_type,amijo,filter] = check_suffDecrease(obj,alpha,x_star,x_k,p0,theta_xk,theta_x_k1, f_x_k1,f_xk,L_k1,L_k,dual_k,filter);
        SuffDecreaseTimeMeas =  toc(measSuffDecreaseTime);

        % either sufficient decrease in cost or constraint violation
        if suffDecrease_flag
            % accept and leave while loop
            feas_res_flag = 0;
            break
        end

        % invoke soc to avoid maratos effect, if necessary
        if ~suffDecrease_flag && alpha == 1  && theta_x_k1 >= theta_xk && obj.opts.enable_SOC

            measSocTime = tic;
            % compute corrected search direction, compute new trial point and check for filter acceptance
            [x_k1,suffConSoC,f_type,amijo] = second_order_correction(obj,x_k,x_star,p0,Bk,args,filter,alpha,theta_xk,f_xk,dkl,L_k1,L_k,dual_k);
            SocTimeMeas = toc(measSocTime);
            % leave while loop if corrected step is acceptable to filter
            if suffConSoC
                feas_res_flag = 0;
                break % means we found a solution with SOC; leave loop
            end

        end

    else % not acceptable to filter

        % invoke soc to avoid maratos effect, if necessary
        if ~filter_Acceptance && alpha == 1 && theta_x_k1 >= theta_xk && obj.opts.enable_SOC

            measSocTime = tic;
            % compute corrected search direction, compute new trial point and check for filter acceptance
            [x_k1,suffConSoC,f_type,amijo] = second_order_correction(obj,x_k,x_star,p0,Bk,args,filter,alpha,theta_xk,f_xk,dkl,L_k1,L_k,dual_k);
            SocTimeMeas = toc(measSocTime);
            if suffConSoC
                feas_res_flag = 0;
                break % means we found a solution with SOC; leave loop
            end

        end

    end

    % if not acceptable to filter and/or no sufficient progress update alpha
    alpha = 1/2*alpha;

    % compute alpha_min
    alpha_min = compute_alpha_min(obj,x_k,dk,p0,theta_xk,filter);

    % if alpha < alpha_min  % just avoid very small step lengths
    %     % invoke feasibility restoration if step-length is below minimum
    %     feas_res_flag = 2;
    %     break
    % end

    if alpha_min == 0 && ~filter_Acceptance
        % invoke feasibility restoration if step-length is below minimum
        feas_res_flag = 2;
        break
    end

end % end of backtracking linesearch
totalBackstepTimeMeas = toc(measBackStepTime);


% only augment if either f-type or amijo are not fulfilled
if (~f_type || ~amijo)

    if LangrangeFilter
        % augment filter with a small margin
        filter = [filter;[theta_xk*(1-gamma_theta), L_k - gamma_phi*theta_xk]];
    else
        % augment filter with a small margin
        filter = [filter;[theta_xk*(1-gamma_theta), f_xk - gamma_phi*theta_xk]];
    end

end

% output of current iterate

sol_iter.x_k1       = x_k1;
sol_iter.dual_k1    = dual_k1;
sol_iter.theta_x_k1 = theta_x_k1;
sol_iter.f_x_k1     = f_x_k1;
sol_iter.alpha_k    = alpha;
sol_iter.dual_qp    = dual_star;
info{iter}.constraint_violation = sol_convio;

measHessApp = tic;
% update quasi-Newton/exact Hessian
if strcmpi(obj.opts.hessian_approx,'BFGS')
    % damped BFGS
    Bk = damped_BFGS(obj,Bk,x_k,p0,sol_iter,iter);

elseif strcmpi(obj.opts.hessian_approx,'gerschgorin')
    % see Betts book
    H  = full(obj.eval_Hessian(x_k1,p0,dual_k1));
    Bk = casos.package.solvers.SequentialCommon.regularize_gerschgorinBound(H);

elseif strcmpi(obj.opts.hessian_approx,'scaledFrob')
    % scaled Frobenius
    H = full(obj.eval_Hessian(x_k1,p0,dual_k1));
    Bk = casos.package.solvers.SequentialCommon.regularize_scaledFrobNorm(H);
elseif strcmpi(obj.opts.hessian_approx,'minFrob')
    % min. Frobenius
    H = full(obj.eval_Hessian(x_k1,p0,dual_k1));
    Bk = casos.package.solvers.SequentialCommon.regularize_minFrobNorm(H);

elseif strcmpi(obj.opts.hessian_approx,'mirroring')
    % Verschueren
    H = full(obj.eval_Hessian(x_k1,p0,dual_k1));
    [L,D] = ldl(H);

    D(D<0) = 1e-6;
    Bk = L*abs(D)*L';

else
    error('Unknown option "%s" for Hessian approximation.',obj.opts.hessian_approx)
end

HessApproxTimeMeas = toc(measHessApp);

% get timing stats
info{iter}.timeStats.HessApproxTime         = HessApproxTimeMeas;
info{iter}.timeStats.totalBackstepTime      = totalBackstepTimeMeas;
info{iter}.timeStats.FilterAcceptTime       = FilterAcceptTimeMeas;
info{iter}.timeStats.SuffDecreaseTime       = SuffDecreaseTimeMeas;

if ~isempty(SocTimeMeas)
    info{iter}.timeStats.SocTime                = SocTimeMeas;
else
    info{iter}.timeStats.SocTime                = 0;
end


end
