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

A = data.A(idx,:);
b = data.b(idx);

% add cones for the clarabel wrapper
if K.z > 0
cones = zeroConeT(K.z);
else
    cones = [];
end



for jj = 1:length(K.s)

    cones = [cones;PSDTriangleConeT(K.s(jj))];

end

% call clarabel mex
sol = clarabel_mex(P,c,A,b,cones,opts);

x  = sol.x;
z_ = sol.z;
s_ = sol.s;


% assign solution
z = sparse(idx,1,z_,m,1);
s = sparse(idx,1,s_,m,1);

if strcmp(sol.status,'Solved')
        % success
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
else
        % primal unbounded / dual infeasible
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    assert(~obj.opts.error_on_fail,'Conic problem is dual infeasible.')
end


    % 
    % // solution status (enumeration) to string for output
    % const char *status_str;
    % switch (solution.status) {
    %     case SolverStatus::Unsolved:
    %         status_str = "Unsolved";
    %         break;
    %     case SolverStatus::Solved:
    %         status_str = "Solved";
    %         break;
    %     case SolverStatus::PrimalInfeasible:
    %         status_str = "PrimalInfeasible";
    %         break;
    %     case SolverStatus::DualInfeasible:
    %         status_str = "DualInfeasible";
    %         break;
    %     case SolverStatus::AlmostSolved:
    %         status_str = "AlmostSolved";
    %         break;
    %     case SolverStatus::AlmostPrimalInfeasible:
    %         status_str = "AlmostPrimalInfeasible";
    %         break;
    %     case SolverStatus::AlmostDualInfeasible:
    %         status_str = "AlmostDualInfeasible";
    %         break;
    %     case SolverStatus::MaxIterations:
    %         status_str = "MaxIterations";
    %         break;
    %     case SolverStatus::MaxTime:
    %         status_str = "MaxTime";
    %         break;
    %     case SolverStatus::NumericalError:
    %         status_str = "NumericalError";
    %         break;
    %     case SolverStatus::InsufficientProgress:
    %         status_str = "InsufficientProgress";
    %         break;
    %     default:
    %         status_str = "Unknown";
    %         break;
    % }
    % 


%     % primal unbounded / dual infeasible
%     obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
%     assert(~obj.opts.error_on_fail,'Conic problem is dual infeasible.')
% elseif status_val == -2
%     % primal infeasible / dual unbounded
%     obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
%     assert(~obj.opts.error_on_fail,'Conic problem is primal infeasible.')
% elseif ismember(status_val, [2 -6 -7])
%     % inaccurate solution
%     obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_NAN;
%     assert(~obj.opts.error_on_fail,'Optimizer did not reach desired accuracy (Status: %s).', obj.info.status)
% elseif status_val < -2
%     % failure
%     obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
%     assert(~obj.opts.error_on_fail,'Optimizer failed (Status: %s).', obj.info.status)
% else
%     % success
%     obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
% end

% parse solution
argout = call(obj.ghan,[argin {x z s}]);

end
