function buildproblem(obj,prob,data,opts)
% Convert SDP structure into conic problem description.
%
% The high-level SDP interface has the form
%
%   min f(x,p) s.t. g(x,p) in Kc, x in Kx
%
% where f = 1/2 x'*H(p)*x + g(p)'*x + d(p) and g = A(p)*x + b(p).
%
% The low-level conic interface has the form
%
%   min 1/2 x'*H*x + g'*x s.t. A*x in Kc, x in Kx
%

sz_x = size_in(prob,0);
sz_p = size_in(prob,1);

% decision variables & parameters
x = casadi.MX.sym('x',sz_x);
p = casadi.MX.sym('p',sz_p);

% sdp problem
[fx,gx] = prob(x,p);

% conic problem data
[Hp,gp,Ap,bp] = data(0,p);

sz_g = size(gx);

% ensure cone has default value
if ~isfield(opts,'Kx')
    opts.Kx.lin = prod(sz_x);
elseif ~isfield(opts.Kx,'lin')
    opts.Kx.lin = 0; 
end
if ~isfield(opts,'Kc')
    opts.Kc.lin = prod(sz_g);
elseif ~isfield(opts.Kc,'lin')
    opts.Kc.lin = 0;
end

% linear variables
Nx_l = opts.Kx.lin;
% cone variables
Nx_c = prod(sz_x) - Nx_l;

% linear constraints
Ng_l = opts.Kc.lin;
% cone constraints
Ng_c = prod(sz_g) - Ng_l;

% separate constraints
Bp = mat2cell(bp,[Ng_l Ng_c],1);

% embed low-level conic solver
% (h,g,a,lba,uba,cba,lbx,ubx,cbx,x0,lam_x0,lam_a0)->(x,cost,lam_a,lam_x)
%
% into high-level SDP interface
% (x0,p,lbx,ubx,cbx,lbg,ubg,cbg,lam_x0,lam_g0)->(x,f,g,lam_x,lam_g,lam_p)
x0 = casadi.MX.sym('x0',sz_x);
lbx = casadi.MX.sym('lbx',Nx_l,1);
ubx = casadi.MX.sym('ubx',Nx_l,1);
cbx = casadi.MX.sym('cbx',Nx_c,1);
lbg = casadi.MX.sym('lbg',Ng_l,1);
ubg = casadi.MX.sym('ubg',Ng_l,1);
cbg = casadi.MX.sym('cbg',Ng_c,1);
lam_x0 = casadi.MX.sym('lam_x0',sz_x);
lam_g0 = casadi.MX.sym('lam_g0',sz_g);

cost = casadi.MX.sym('cost');
lam_a = casadi.MX.sym('lam_a',sz_g);
lam_x = casadi.MX.sym('lam_x',sz_x);

% options for Casadi functions
% NOTE: As per Casadi v3.6.5, functions' input and output names 
% must be mutually exclusive unless explicity permitted.
fopt = struct('allow_duplicate_io_names',true);

% input function
obj.fhan = casadi.Function('f', ...
            {x0 p lbx ubx cbx lbg ubg cbg lam_x0 lam_g0}, ...
            {Hp gp Ap Bp{1}+lbg Bp{1}+ubg Bp{2}+cbg lbx ubx cbx x0 lam_x0 lam_g0}, ...
            {'x0' 'p' 'lbx' 'ubx' 'cbx' 'lbg' 'ubg' 'cbg' 'lam_x0' 'lam_g0'}, ...
            {'h' 'g' 'a' 'lba' 'uba' 'cba' 'lbx' 'ubx' 'cbx' 'x0' 'lam_x0' 'lam_a0'}, ...
            fopt ...
);

% output function
obj.ghan = casadi.Function('g', ...
            {p x cost lam_a lam_x}, ...
            {x fx gx lam_x lam_a nan(size(p))}, ...
            {'p' 'x' 'cost' 'lam_a' 'lam_x'}, ...
            {'x' 'f' 'g' 'lam_x' 'lam_g' 'lam_p'}, ...
            fopt ...
);

% number of decision variables
obj.sdp_info.nx = numel(sz_x);

% number of constraints (lin/conic)
obj.sdp_info.nc.lin = Ng_l;
obj.sdp_info.nc.conic = Ng_c;

% number of nonzero elements in linear constraints (Ap)
obj.sdp_info.nnz_lin = sparsity(Ap).nnz;


% TODO:
% Embed conic interface into (expr_in)->(expr_out).
% f = obj.fhan;
% g = obj.ghan;
% 
% % low-level interface: g.S.f
% % high-level interface: out.(g.S.f).in = (out.g).S.(f.in)
% obj.fhan = casadi.Function(name(f),mx_in(fun_in),call(f,mx_out(fun_in)),name_in(fun_in),name_out(f));
% obj.ghan = casadi.Function(name(g),mx_in(g),call(fun_out,mx_out(g)),name_in(g),name_out(fun_out));

end
