function argout = eval(obj,argin)
% Call SCS interface.

prob = call(obj.fhan,cell2struct(argin',fieldnames(obj.args_in)));
cone = call(obj.cone,cell2struct(argin',fieldnames(obj.args_in)));

% to double
data.P = sparse(prob.P);
data.A = sparse(prob.A);
data.b =  full(prob.b);
data.c =  full(prob.c);
% cone
K = structfun(@full,cone,'UniformOutput',false);

% options to SCS
opts = obj.sdpopt.solveroptions;

% call SCS
[x,y,s,info] = scs(data,K,opts);

if ~obj.sdpopt.error_on_fail
    % continue regardless of feasibility
elseif info.status_val == -1
    % primal unbounded / dual infeasible
    error('Conic problem is dual infeasible.')
elseif info.status_val == -2
    % primal infeasible / dual unbounded
    error('Conic problem is primal infeasible.')
elseif ismember(info.status_val, [2 -6 -7])
    % inaccurate solution
    error('Optimizer did not reach desired accuracy (Status: %s).', info.status)
elseif info.status_val < -2
    % failure
    error('Optimizer failed (Status: %s).', info.status)
end

% parse solution
argout = call(obj.ghan,[argin {x y s}]);

end
