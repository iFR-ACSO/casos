function argout = eval(obj,argin)
% Call Clarabel interface.

% evaluate problem structure
prob = call(obj.fhan,cell2struct(argin',fieldnames(obj.args_in)));
cone = call(obj.cone,cell2struct(argin',fieldnames(obj.args_in)));

% to double
P = sparse(prob.P);
data.A = sparse(prob.A);
data.b =  full(prob.b);
c =  full(prob.c);
% cone
K = structfun(@full,cone,'UniformOutput',false);

% options to SCS
opts = obj.opts.Clarabel;
% disable verbosity by default
% if ~isfield(opts,'verbose'), opts.verbose = 0; end

% joint lower bounds
% -A(x) + s = 0, s in [lb ub]
lb = sparse(K.bl);
ub = sparse(K.bu);

ml = length(lb)+1;
m  = length(data.b);

% reoder slack variables
%idx = 1:length(data.b);

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

A = data.A(idx,:);
b = data.b(idx);

% Initialize an empty cell array for cones
cones = cell(1, 1); 

% Add zero cones for the Clarabel wrapper
if K.z > 0
    cones{end+1} = zeroConeT(K.z);
end

% Add second order cones
if K.q > 0
    cones{end+1} = SecondOrderConeT(K.q);
end

% Stack PSD cones
num_psd = length(K.s);
if num_psd > 0
    psd_cones = arrayfun(@PSDTriangleConeT, K.s, 'UniformOutput', false);
    % concatenate
    cones = [cones, psd_cones']; 
end

% Convert cell array to a proper structure 
cones = vertcat(cones{:});

% call clarabel mex
sol = clarabel_mex(P,c,A,b,cones,opts);

% extract primal, dual and slack variables
x  = sol.x;
z_ = sol.z;
s_ = sol.s;

% assign solution
z = sparse(idx,1,z_,m,1);
s = sparse(idx,1,s_,m,1);

% get solution status
if strcmp(sol.status,'Solved')
    % success
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
elseif strcmp(sol.status,'AlmostSolved')
    % success because we are feasible but did not solve until thresholds
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
elseif strcmp(sol.status,'Unsolved')
    % unsolved
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
elseif strcmp(sol.status,'MaxIterations')
    % max iterations reached
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
elseif strcmp(sol.status,'MaxTime')
    % max time reached
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
elseif strcmp(sol.status,'InsufficientProgress')
    % InsufficientProgress
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
elseif strcmp(sol.status,'Unknown')
    % InsufficientProgress
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
elseif strcmp(sol.status,'NumericalError')
    % InsufficientProgress
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
else
    % primal unbounded / dual infeasible
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    assert(~obj.opts.error_on_fail,['Conic problem is ' sol.status])
end

% parse solution
argout = call(obj.ghan,[argin {x z s}]);

end
