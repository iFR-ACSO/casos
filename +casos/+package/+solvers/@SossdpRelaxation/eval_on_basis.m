function argout = eval_on_basis(obj,argin)
% Evaluate SDP relaxation to solve SOS problem.

% project arguments to obtain SDP inputs
% only linear coefficients are handled (p, lbx, ubx, lbg, ubg)
% TODO: handle SOS inputs
[~,p,lbx,ubx,~,lbg,ubg] = argin{:};

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

% store Gram representations
obj.info.gram = call(obj.sdp2gram,struct('x_sol',sdpsol.x));


% build polynomial solution
argout = {sossol.x sossol.f sossol.g sossol.lam_x sossol.lam_g};

end
