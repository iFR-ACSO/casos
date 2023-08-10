function buildproblem(obj,prob,data)
% Convert SDP structure into conic problem description.
%
% The high-level SDP interface has the form
%
%   min f(x,p) s.t. g(x) in Kg, x in Kx
%
% where f = 1/2 x'*H(p)*x + g(p)'*x + d(p) and g = A(p)*x + b(p).
%
% The low-level conic interface has the form
%
%   min 1/2 x'*H*x + g'*x s.t. A*x in Ka, x in Kx
%

sz_x = size_in(prob,0);
sz_p = size_in(prob,1);

% decision variables & parameters
x = casadi.MX.sym('x',sz_x);
p = casadi.MX.sym('p',sz_p);

% sdp problem
[fx,gx] = prob(x,p);

% conic problem data
[Hp,gp,Ap,bp] = data(p);

sz_g = size(gx);

% embed low-level conic solver
% (h,g,a,lba,uba,lbx,ubx,x0,lam_x0,lam_a0)->(x,cost,lam_a,lam_x)
%
% into high-level SDP interface
% (x0,p,lbx,ubx,lbg,ubg,lam_x0,lam_g0)->(x,f,g,lam_x,lam_g,lam_p)
x0 = casadi.MX.sym('x0',sz_x);
lbx = casadi.MX.sym('lbx',sz_x);
ubx = casadi.MX.sym('ubx',sz_x);
lbg = casadi.MX.sym('lbg',sz_g);
ubg = casadi.MX.sym('ubg',sz_g);
lam_x0 = casadi.MX.sym('lam_x0',sz_x);
lam_g0 = casadi.MX.sym('lam_g0',sz_g);

cost = casadi.MX.sym('cost');
lam_a = casadi.MX.sym('lam_a',sz_g);
lam_x = casadi.MX.sym('lam_x',sz_x);

% input function
obj.fhan = casadi.Function('f', ...
            {x0 p lbx ubx lbg ubg lam_x0 lam_g0}, ...
            {Hp gp Ap bp+lbg bp+ubg lbx ubx x0 lam_x0 lam_g0}, ...
            {'x0' 'p' 'lbx' 'ubx' 'lbg' 'ubg' 'lam_x0' 'lam_g0'}, ...
            {'h' 'g' 'a' 'lba' 'uba' 'lbx' 'ubx' 'x0' 'lam_x0' 'lam_a0'} ...
);

% output function
obj.ghan = casadi.Function('g', ...
            {x cost lam_a lam_x}, ...
            {x fx gx lam_x lam_a nan(size(p))}, ...
            {'x' 'cost' 'lam_a' 'lam_x'}, ...
            {'x' 'f' 'g' 'lam_x' 'lam_a' 'lam_p'} ...
);

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
