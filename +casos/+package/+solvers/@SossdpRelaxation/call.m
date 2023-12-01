function argout = call(obj,argin)
% Call SDP relaxation to solve SOS problem.

% select Gram decision variables
Ix = [false(size(obj.monom_xl,1),1); true(size(obj.gram_x,1),1); false(size(obj.gram_g,1),1)];
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
sol_xl = sdpsol.x(find(~Ix & ~Ig));
sol_xs = sdpsol.x(find(Ix));
sol_gs = sdpsol.x(find(Ig));
sol_gl = sdpsol.g(find(~Jg));
sol_lam_xl = sdpsol.lam_x(find(~Ix & ~Ig));
sol_lam_xs = sdpsol.lam_x(find(Ix));
sol_lam_gs = sdpsol.lam_x(find(Ig));
sol_lam_gl = sdpsol.lam_g(find(~Jg));

% build polynomial solution
argout{2} = casos.PS(obj.monom_f,sol_f); % f
argout{1} = [casos.PS(obj.monom_xl,sol_xl); casos.PS(obj.gram_x,sol_xs)]; % x
argout{3} = [casos.PS(obj.monom_gl,sol_gl); casos.PS(obj.gram_g,sol_gs)]; % g
argout{4} = [casos.PS(obj.monom_xl,sol_lam_xl); casos.PS(obj.gram_x,sol_lam_xs)]; % lam_x
argout{5} = [casos.PS(obj.monom_gl,sol_lam_gl); casos.PS(obj.gram_g,sol_lam_gs)]; % lam_g

end
