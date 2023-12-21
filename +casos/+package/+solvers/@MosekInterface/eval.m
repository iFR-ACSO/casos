function argout = eval(obj,argin)
% Call MOSEK interface.

msk_prob = obj.cone;

% evaluate problem structure
prob = call(obj.fhan,cell2struct(argin',fieldnames(obj.args_in)));
data = call(obj.barv,cell2struct(argin',fieldnames(obj.args_in)));

% to double
msk_prob.a = sparse(prob.a);
msk_prob.c = sparse(prob.c);
msk_prob.f = sparse(prob.f);
msk_prob.g = full(prob.g);
msk_prob.blc = full(prob.blc);
msk_prob.buc = full(prob.buc);
msk_prob.blx = full(prob.blx);
msk_prob.bux = full(prob.bux);
% structures
msk_prob.bara.val = full(data.a);
msk_prob.barc.val = full(data.c);
msk_prob.barf.val = full(data.f);
%TODO: initial guess, quadratic cost

% options to MOSEK
msk_param = obj.opts.mosek_param;
msk_echo  = obj.opts.mosek_echo;

% build MOSEK command string
msk_cmd = sprintf('minimize echo(%d) info statuskeys(0)',msk_echo);

% call MOSEK
[rcode,res] = mosekopt(msk_cmd,msk_prob,msk_param);

% store info (if any)
if isfield(res,'info')
    obj.info = res.info;
end

% pre-initialize solution struct
sol = struct('pobjval',0,'xx',0,'barx',0,'slc',0,'suc',0,'slx',0,'sux',0,'doty',0,'bars',0);

% retrieve solution
if isfield(res,'sol') && isfield(res.sol,'itr')
    % interior-point solution available
    msk_sol = res.sol.itr;
elseif isfield(res,'sol') && isfield(res.sol,'bas')
    % solution from simplex available
    msk_sol = res.sol.bas;
else
    % no solution available (some error occured)
    msk_sol = [];
end

% check solution
if ~isempty(msk_sol)
    % check problem status
    msk_prosta = msk_sol.prosta;
    switch (msk_prosta)
        case {'PRIMAL_AND_DUAL_FEASIBLE' 'PRIMAL_FEASIBLE' 'DUAL_FEASIBLE'}
            % feasible problem
            msk_solsta = msk_sol.solsta;
            % check solution status
            if ismember(msk_solsta,{'OPTIMAL' msk_prosta})
                % solution status matches 
                obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
            else
                % solution status infeasible or unknown
                obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
                assert(~obj.opts.error_on_fail,'Problem appears feasible (%s) but solution is not (%s).',msk_prosta,msk_solsta)
            end
        case {'PRIMAL_INFEASIBLE' 'DUAL_INFEASIBLE' 'PRIMAL_AND_DUAL_INFEASIBLE' 'PRIMAL_INFEASIBLE_OR_UNBOUNDED'}
            % infeasible problem
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
            assert(~obj.opts.error_on_fail,'Problem is primal and/or dual infeasible or unbounded.')
        otherwise
            % problem status unknown or ill-posed
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
            assert(~obj.opts.error_on_fail,'Problem status unknown (might be ill-posed).')
    end
    
    % retrieve solution
    sol.pobjval = msk_sol.pobjval;
    sol.xx = msk_sol.xx;
    sol.barx = msk_sol.barx;
    sol.slc = msk_sol.slc;
    sol.suc = msk_sol.suc;
    sol.slx = msk_sol.slx;
    sol.sux = msk_sol.sux;
    sol.doty = msk_sol.doty;
    sol.bars = msk_sol.bars;

% check return code
elseif rcode > 0
    % unexpected termination, warning, or error
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
    assert(~obj.opts.error_on_fail,'Call to mosekopt failed: %s (%s).',res.rmsg,res.rcodestr);
else
    % unexpected error (no solution)
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
    assert(~obj.opts.error_on_fail,'Call to mosekopt did not yield a solution.')
end

% parse solution
argout = call(obj.ghan,struct2cell(sol));

end
