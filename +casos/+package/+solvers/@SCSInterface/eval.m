function argout = eval(obj,argin)
% Call SCS interface.

% evaluate problem structure
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

% joint lower bounds
% -A(x) + s = 0, s in [lb ub]
lb = sparse(K.bl);
ub = sparse(K.bu);

ml = length(lb)+1;
m  = length(data.b);

% reoder slack variables
idx = 1:length(data.b);

% remove trivial constraints
J = false(size(data.b));

% determine equality constraints
I0 = find(lb == ub);
% remove box constraint
J(1+I0) = true;
% set equality constraint
% -A(x) + s = -b, s = 0
data.b(1+I0) = -lb(I0);

% determine right-open constraints
Ipos = find(isinf(ub) & ~isinf(lb));
% remove box constraint
J(1+Ipos) = true;
% set lower bound
% -A(x) + s = -lb, s >= 0  <->  A(x) = lb + s
data.b(1+Ipos) = -lb(Ipos);

% determine left-open constraints
Ineg = find(isinf(lb) & ~isinf(ub));
% remove box constraint
J(1+Ineg) = true;
% invert constraint
data.A(1+Ineg) = -data.A(1+Ineg);
% set upper bound
% +A(x) + s = ub, s >= 0  <->  A(x) = ub - s
data.b(1+Ineg) = ub(Ineg);

% determine unbounded constraints
Iinf = find(isinf(lb) & isinf(ub));
% remove box constraint
J(1+Iinf) = true;

% modify cone
K.z = length(I0);
K.l = length(Ipos) + length(Ineg);
K.bl(J(2:ml)) = [];
K.bu(J(2:ml)) = [];

% remove box constraint slack variable
J(1) = all(J(2:ml));

% reorder constraints
idx = [1+[I0; Ipos; Ineg]; find(~J)];

data.A = data.A(idx,:);
data.b = data.b(idx);

% call SCS
[x,y_,s_,obj.info] = scs(data,K,opts);

% assign solution
y = sparse(idx,1,y_,m,1);
s = sparse(idx,1,s_,m,1);

status_val = obj.info.status_val;
if status_val == -1
    % primal unbounded / dual infeasible
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    assert(~obj.sdpopt.error_on_fail,'Conic problem is dual infeasible.')
elseif status_val == -2
    % primal infeasible / dual unbounded
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    assert(~obj.sdpopt.error_on_fail,'Conic problem is primal infeasible.')
elseif ismember(status_val, [2 -6 -7])
    % inaccurate solution
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_NAN;
    assert(~obj.sdpopt.error_on_fail,'Optimizer did not reach desired accuracy (Status: %s).', obj.info.status)
elseif status_val < -2
    % failure
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
    assert(~obj.sdpopt.error_on_fail,'Optimizer failed (Status: %s).', obj.info.status)
else
    % success
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
end

% parse solution
argout = call(obj.ghan,[argin {x y s}]);

end
