function argout = eval_on_basis(obj,argin)
% Evaluate LMI to solve SOS problem.

% project arguments to obtain LMI parameters
% only linear coefficients are handled (p, lbx, ubx, lbg, ubg)
% TODO: handle SOS inputs
[~,p,lbx,ubx,~,lbg,ubg] = argin{:};

% evaluate problem structure
prob = call(obj.fhan,{p lbx ubx lbg ubg});
% retrieve cones
cone = obj.cone;

% to double
A = sparse(prob{1});
b = sparse(prob{2});
c = sparse(prob{3});
lbx = sparse(lbx);
ubx = sparse(ubx);
lbg = sparse(lbg);
ubg = sparse(ubg);

% dimensions of original problem
nl = length(lbx);
ml = length(lbg);

% detect infinite bounds
Iinf = isinf([lbx; ubx; lbg; ubg]);

% remove trivial constraints from LMI
I = [Iinf; false(length(c)-2*nl-2*ml,1)];

% detect equality constraints
Ilx = find(lbx == ubx);
Ilg = find(lbg == ubg);
% remove lower bounds
I(Ilx) = true;
I(2*nl+Ilg) = true;
% remove upper bounds
I(nl+Ilx) = true;
I(2*nl+ml+Ilg) = true;

% reorder variables in LMI
idx = [find(~I) Ilx 2*nl+Ilg];

% reorder and purge problem structure
A = A(:,idx);
c = c(idx);

% modify cone
cone{1}.dim = cone{1}.dim - nnz(I);
cone{end+1} = struct('type','free','dim',length(Ilx)+length(Ilg));

% options to alfonso
opts = obj.opts.alfonso;
% disable output by default
if ~isfield(opts,'verbose'), opts.verbose = 0; end

% call alfonso
res = alfonso_simple(c,A,b,cone,[],opts);

% assign full solution
x = sparse(idx,1,res.x,length(I),1);
s = sparse(idx,1,res.s,length(I),1);
% copy slack variables for infinite lower bounds
Iinf_lbg = find(isinf(lbg));
s(2*nl+Iinf_lbg) = -s(2*nl+ml+Iinf_lbg);
lbg(Iinf_lbg) = -ubg(Iinf_lbg);

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
argout = call(obj.ghan,{res.dObj x s res.y lbg});

end
