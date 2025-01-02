function argout = eval_on_basis(obj,argin)
% Evaluate LMI to solve SOS problem.

% project arguments to obtain LMI parameters
% only linear coefficients are handled (p, lbx, ubx, lbg, ubg)
% TODO: handle SOS inputs
[~,p,lbx,ubx,~,lbg,ubg] = argin{:};

% evaluate problem structure
prob = call(obj.fhan,struct('p',p));

% to double
A = sparse(prob.A);
b = sparse(prob.b);
c = sparse(prob.c);

% options to alfonso
opts = obj.opts.alfonso;
% disable output by default
if ~isfield(opts,'verbose'), opts.verbose = 0; end

% call alfonso
res = alfonso_simple(c,A,b,obj.cone,[],opts);

if res.status > 0
    % success
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
elseif any(res.status == [-2 -3])
    % LMI is dual infeasible
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    assert(~obj.opts.error_on_fail,'Sum-of-squares problem is primal infeasible.')
elseif res.status == -1
    % LMI is primal infeasible
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
    assert(~obj.opts.error_on_fail,'Sum-of-squares problem is dual infeasible.')
else
    % status unknown (might be ill-posed)
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
    assert(~obj.opts.error_on_fail,'Problem status unknown (might be ill-posed).')
end

% store alfonso infos (TODO)
obj.info.alfonso_status = res.statusString;

% build polynomial solution
argout = call(obj.ghan,{res.dObj res.x res.s res.y});

end
