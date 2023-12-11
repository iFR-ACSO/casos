function argout = call(obj,argin)
% Call SDP relaxation to solve SOS problem.

% select Gram decision variables
Ig = [false(size(obj.monom_xl,1),1); false(size(obj.gram_x,1),1); true(size(obj.gram_g,1),1)];
% select SOS constraints
Jg = [false(size(obj.monom_gl,1),1); true(size(obj.monom_gs,1),1)];

% project arguments to obtain SDP inputs
% only linear coefficients are handled (p, lbx, ubx, lbg, ubg)
% TODO: handle SOS inputs
p = poly2basis(argin{2}, obj.monom_p);
lbx = poly2basis(argin{3}, obj.monom_xl);
ubx = poly2basis(argin{4}, obj.monom_xl);
lbg = poly2basis(argin{6}, obj.monom_gl);
ubg = poly2basis(argin{7}, obj.monom_gl);

% prepare arguments to SDP
args.p = casadi.DM(p);
args.lbx = casadi.DM(lbx);
args.ubx = casadi.DM(ubx);
args.lbg = casadi.DM([lbg; zeros(nnz(Jg),1)]);
args.ubg = casadi.DM([ubg; zeros(nnz(Jg),1)]);

% call SDP solver
sdpsol = call(obj.sdpsolver, args);

% select solutions
sol_f  = sdpsol.f;
sol_x  = sdpsol.x(find(~Ig));
sol_g  = [
    sdpsol.x(find(Ig))
    sdpsol.g(find(~Jg))
];
sol_lam_x = sdpsol.lam_x(find(~Ig));
sol_lam_g = [
    sdpsol.lam_x(find(Ig))
    sdpsol.lam_g(find(~Jg))
];

% build polynomial solution
argout{2} = casos.PS(obj.monom_f,sol_f); % f
argout{1} = casos.PS(obj.basis_x_out,sol_x); % x
argout{3} = casos.PS(obj.basis_g_out,sol_g); % g
argout{4} = casos.PS(obj.basis_x_out,sol_lam_x); % lam_x
argout{5} = casos.PS(obj.basis_g_out,sol_lam_g); % lam_g

end
