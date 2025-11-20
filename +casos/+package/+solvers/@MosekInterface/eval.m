function argout = eval(obj,argin)
% Call MOSEK interface.

msk_prob = obj.cone;

args = cell2struct(argin',fieldnames(obj.args_in));

% compute cholseky numerically
if nnz(obj.args_in.h) > 0 && strcmp(obj.opts.cholesky_method,'numerical')
    args.h = chol((args.h));
end


% evaluate problem structure
prob = call(obj.fhan,args);
data = call(obj.barv,args);

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
msk_param.MSK_IPAR_AUTO_UPDATE_SOL_INFO = 'MSK_ON';
msk_param.MSK_IPAR_INTPNT_BASIS         ='MSK_BI_NEVER';
% msk_param.MSK_IPAR_INTPNT_STARTING_POINT = 'MSK_STARTING_POINT_CONSTANT';
msk_echo  = obj.opts.mosek_echo;
    

% build MOSEK command string
msk_cmd = sprintf('minimize echo(%d) info statuskeys(0)',msk_echo);

% call MOSEK
[rcode,res] = mosekopt(msk_cmd,msk_prob,msk_param);

% store info (if any)
if isfield(res,'info')
    obj.info.mosek_info = res.info;
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
    obj.info.mosek_prosta = msk_sol.prosta;
    obj.info.mosek_solsta = msk_sol.solsta;
    switch (msk_sol.prosta)
        case {'PRIMAL_AND_DUAL_FEASIBLE' 'PRIMAL_FEASIBLE' 'DUAL_FEASIBLE'}
            % feasible problem, check solution status
            if ismember(msk_sol.solsta,{'OPTIMAL' msk_sol.prosta})
                % solution status matches 
                obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
            else
                % solution status infeasible or unknown
                obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
                assert(~obj.opts.error_on_fail,'Problem appears feasible (%s) but solution is not (%s).',msk_sol.prosta,msk_sol.solsta)
            end
        case {'PRIMAL_INFEASIBLE' 'DUAL_INFEASIBLE' 'PRIMAL_AND_DUAL_INFEASIBLE' 'PRIMAL_INFEASIBLE_OR_UNBOUNDED'}
            % infeasible problem
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
            assert(~obj.opts.error_on_fail,'Problem is primal and/or dual infeasible or unbounded.')
         case {'UNKNOWN'}
            % in this case mosek can not decide about the solution so we
            % need some additional check here!

             % Compute relative primal-dual gap
            relativeGap = abs(res.info.MSK_DINF_SOL_ITR_PRIMAL_OBJ - res.info.MSK_DINF_SOL_ITR_DUAL_OBJ) ...
                          / max([abs(res.info.MSK_DINF_SOL_ITR_PRIMAL_OBJ), abs(res.info.MSK_DINF_SOL_ITR_DUAL_OBJ)]);
            
            % Compute max primal constraint violation across different categories
            primMaxVio = max([abs(res.info.MSK_DINF_SOL_ITR_PVIOLCON), ...   % General constraints
                              abs(res.info.MSK_DINF_SOL_ITR_PVIOLACC), ...   % Affine conic constraints
                              abs(res.info.MSK_DINF_SOL_ITR_PVIOLCONES), ... % Conic constraints
                              abs(res.info.MSK_DINF_SOL_ITR_PVIOLVAR)]);     % Variable bound violations
            
            % Compute max dual constraint violation across different categories
            dualMaxVio = max([abs(res.info.MSK_DINF_SOL_ITR_DVIOLCON), ... % General constraints
                              abs(res.info.MSK_DINF_SOL_ITR_DVIOLACC), ... % Affine conic constraints
                              abs(res.info.MSK_DINF_SOL_ITR_DVIOLCONES), ... % Conic constraints
                              abs(res.info.MSK_DINF_SOL_ITR_DVIOLVAR)]); % Variable bound violations
            
            % Compute optimality measure from Mosek's solver status
            % It should converge to +1 to be optimal
            optMeas = res.info.MSK_DINF_INTPNT_OPT_STATUS > obj.opts.mosek_augmented_feasCheck.optMeas;
            
            % Check for potential ill-conditioning by looking at large norm values
            maxNorm = max([res.info.MSK_DINF_SOL_ITR_NRM_XX, ... % Norm of primal variables
                           res.info.MSK_DINF_SOL_ITR_NRM_Y]);    % Norm of dual variables
            
            % Define an adaptive feasibility tolerance based on problem scale
            % If the objective values are large, allow a slightly relaxed threshold
            feasibilityTol = obj.opts.mosek_augmented_feasCheck.feasTol * max(1, max([abs(res.info.MSK_DINF_SOL_ITR_PRIMAL_OBJ), ...
                                                                                      abs(res.info.MSK_DINF_SOL_ITR_DUAL_OBJ)]));
            
            % log additional decision logic
            obj.info.mosek_acceptable_info.primMaxVio  = primMaxVio;
            obj.info.mosek_acceptable_info.primMaxVio  = dualMaxVio;
            obj.info.mosek_acceptable_info.relativeGap = relativeGap;
            obj.info.mosek_acceptable_info.optMeas     = optMeas;
            obj.info.mosek_acceptable_info.maxNorm     = maxNorm;
            
            % store the decision values 
            obj.info.mosek_acceptable_info.tolerances.primMaxVio  = feasibilityTol;
            obj.info.mosek_acceptable_info.tolerances.dualMaxVio  = feasibilityTol;
            obj.info.mosek_acceptable_info.tolerances.relativeGap = obj.opts.mosek_augmented_feasCheck.relativeGap;
            obj.info.mosek_acceptable_info.tolerances.maxNorm     = obj.opts.mosek_augmented_feasCheck.maxNorm;
            obj.info.mosek_acceptable_info.tolerances.optMeas     = obj.opts.mosek_augmented_feasCheck.optMeas;
            
            % Decision logic for solution acceptability
            if relativeGap < obj.opts.mosek_augmented_feasCheck.relativeGap && ...           % Ensure relative gap is less than 5%
               primMaxVio <= feasibilityTol && ... % Check primal feasibility within adaptive tolerance
               dualMaxVio <= feasibilityTol && ... % Check dual feasibility within adaptive tolerance
               optMeas && ...                       % Ensure optimality measure is reasonable
               maxNorm < obj.opts.mosek_augmented_feasCheck.maxNorm                       % Avoid using solutions with extremely large norms (ill-conditioning)
            
                obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS; % Solution is acceptable
                
            else
                obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED; % Solution is questionable
                assert(~obj.opts.error_on_fail,'Problem status unknown (might be ill-posed).')
            end


            
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
argout = call(obj.ghan,[struct2cell(sol)' argin]);

end
