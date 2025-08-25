%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Short Description:
%
% Check if new iterate is acceptable to filter. Compute constraint
% violation and cost at new iterate. Compare the pair {cost constraint
% vio.} to filter entries. If at least one of them is better in cost or
% constrain violation, than accept.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [x_k1,theta_x_k1,f_x_k1 ,filter_Acceptance,all_violations] = check_filter_acceptance(obj, ...
    filter, ...
    alpha, ...
    x_k, ...
    dual_k, ...
    dk, ...
    dkl, ...
    p0, ...
    args)

% use Langrange function instead of cost function in filter
LangrangeFilter = obj.opts.filter_struct.LangrangeFilter;

% compute new solution candidate
x_k1    = full(x_k     + alpha*dk);
if strcmpi(obj.opts.conVioCheck,'signed-distance')
    % compute constraint violation of new solution candidate
    args_conVio     =  args;
    args_conVio{2}  =  [p0; x_k1];
    args_conVio{3}  = -inf(obj.init_para.conVio.no_con,1);
    args_conVio{4}  =  inf(obj.init_para.conVio.no_con,1);


    sol_convio = eval_on_basis(obj.solver_conVio, args_conVio);


    % extract signed-distances
    all_violations = sol_convio{1};

    % get the largest one
    theta_x_k1 = full(max(0,max(all_violations)));
else
    % evalaute at current solution and samples
    g_val = full(obj.eval_constraintSamples(x_k1,obj.opts.userSample));

    % from all constraints get the smallest function value
    all_violations  = min(g_val);

    % get the overall smallest value
    min_vio = min(min(g_val));
    if min_vio < 0
        theta_x_k1 = abs(min_vio);
    else
        theta_x_k1 = 0; % all samples positive; no violation
    end

end
% cost at trial point
if LangrangeFilter
    dual_k1 = full(dual_k  + alpha*dkl);
    L_k1     = full(obj.eval_L(x_k1,p0,dual_k1));
else
    f_x_k1   = full(obj.eval_cost(x_k1,p0));
end


% check filter acceptance
% get filter inputs
theta_l = filter(:,1);
f_l     = filter(:,2);

% new point lies in forbidden region if both are larger than filter entries
dominance_bool      = [];

% boundary belongs to forbidden region i.e. >=
dominance_bool(:,1) = theta_x_k1 >= theta_l;

if LangrangeFilter
    dominance_bool(:,2) = L_k1     >= f_l;
else
    dominance_bool(:,2) = f_x_k1  >= f_l;
end

% check pairs; if both one, means pair lies in forbidden region
dominance_bool = all(dominance_bool, 2);

if any(dominance_bool) % means not acceptable to filter
    filter_Acceptance = 0;
else
    filter_Acceptance = 1;
end

end