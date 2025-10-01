function buildproblem(obj,nlsos)
% Build additional problem structures for filter-linesearch.

buildproblem@casos.package.solvers.SequentialCommon(obj,nlsos);

%% set up feasibility restoration phase
base_x = sparsity(nlsos.x);

I = true(length(nlsos.g),1);
% for idx = 1:length(nlsos.g)
%     I(idx) = ~is_linear(nlsos.g(idx),nlsos.x);
% end

% get gram half-basis for nonlinear constraints
[~,~,z] = grambasis(nlsos.g,I);

% build unit vectors
base_s0 = gramunit(z);

r  = casos.PS.sym('r',sum(I));
s0 = casos.PD(base_s0);

nlsos_feas.x = [r;nlsos.x];
nlsos_feas.g = nlsos.g(I) + r.*s0;

x_R   = casos.PS.sym('x_R',base_x);

if strcmpi(obj.opts.feasibility_restoration,'simple')
    % we simply mininimize on the constraint manifold
    Phi = [];
    lambda = [];
elseif strcmpi(obj.opts.feasibility_restoration,'regularize')
    % check how far we are from the original problem/solution
    e   = nlsos_feas.x(sum(I)+1:end) - nlsos.x;
    Phi = 1/2*dot(e,e);
    lambda   = casos.PS.sym('l');
elseif   strcmpi(obj.opts.feasibility_restoration,'cost')
    Phi = 1/2*nlsos.f;
    lambda   = casos.PS.sym('l');
else
    error('Unknown option "%s" for feasibility restoration.',obj.opts.feasibility_restoration)
end

nlsos_feas.f = sum(r) + lambda*Phi ;

nlsos_feas.p = [nlsos.p; x_R;lambda];

obj.feasRes_para.n_r = sum(I);
obj.feasRes_para.length_dualOut = length(nlsos_feas.x)-length(nlsos_feas.g);

sosopt.sossol         = obj.opts.sossol;
sosopt.sossol_options = obj.opts.sossol_options;
sosopt.Kx.lin         = length(nlsos_feas.x);
sosopt.Kc.sos         = length(nlsos_feas.g);
sosopt.error_on_fail  = false;
sosopt.verbose        = 1;
sosopt.max_iter       = 100;
obj.feas_res_solver   =  casos.package.solvers.FeasibilityRestoration('feasRes',nlsos_feas,sosopt);

% total build time for both actual problem and feasibility
% restoration
obj.display_para.solver_build_time      = obj.display_para.solver_build_time +obj.feas_res_solver.display_para.solver_build_time;

end
