function argout = call(obj,argin)
% Call SDP relaxation to solve SOS problem.

% project arguments to obtain SDP inputs
% only linear coefficients are handled (p, lbx, ubx, lbg, ubg)
% TODO: handle SOS inputs
p = poly2basis(argin{2}, obj.sparsity_p);
lbx = poly2basis(argin{3}, obj.sparsity_xl);
ubx = poly2basis(argin{4}, obj.sparsity_xl);
lbg = poly2basis(argin{6}, obj.sparsity_gl);
ubg = poly2basis(argin{7}, obj.sparsity_gl);

% prepare arguments to SDP
args.p = casadi.DM(p);
args.lbx = casadi.DM(lbx);
args.ubx = casadi.DM(ubx);
args.lbg = casadi.DM([lbg; sparse(nnz(obj.sparsity_gs),1)]);
args.ubg = casadi.DM([ubg; sparse(nnz(obj.sparsity_gs),1)]);

% call SDP solver
sdpsol = call(obj.sdpsolver, args);

% retrieve SOS coordinates
sossol = call(obj.gram2sos,sdpsol);

% build polynomial solution
argout{1} = casos.package.polynomial(obj.sparsity_x,sossol.x);
argout{2} = casos.package.polynomial(obj.sparsity_f,sossol.f);
argout{3} = casos.package.polynomial(obj.sparsity_g,sossol.g);
argout{4} = casos.package.polynomial(obj.sparsity_x,sossol.lam_x);
argout{5} = casos.package.polynomial(obj.sparsity_g,sossol.lam_g);

end
