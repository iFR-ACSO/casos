function argout = eval_on_basis(obj,argin)
% Evaluate LMI to solve SOS problem.

% project arguments to obtain LMI parameters
% only linear coefficients are handled (p, lbx, ubx, lbg, ubg)
% TODO: handle SOS inputs
[~,p,lbx,ubx,~,lbg,ubg] = argin{:};

% to double
lbx = sparse(lbx);
ubx = sparse(ubx);
lbg = sparse(lbg);
ubg = sparse(ubg);

% detect infinite bounds
Iinf_lbx = isinf(lbx);
Iinf_ubx = isinf(ubx);
Iinf_lbg = isinf(lbg);
Iinf_ubg = isinf(ubg);
Iinf = [Iinf_lbx; Iinf_ubx; Iinf_lbg; Iinf_ubg];

% remove infinite bounds (1)
lbx(Iinf_lbx) = 0;
ubx(Iinf_ubx) = 0;
lbg(Iinf_lbg) = 0;
ubg(Iinf_ubg) = 0;

% prepare arguments to LMI
args = call(obj.lmiargs,struct( ...
                            'p',casadi.DM(p), ...
                            'lbx', casadi.DM(lbx), ...
                            'ubx', casadi.DM(ubx), ...
                            'lbg', casadi.DM(lbg), ...
                            'ubg', casadi.DM(ubg) ...
));

% remove infinite bounds (2)
args.ubx(find(Iinf)) = 0;

% call LMI solver
lmisol = call(obj.sdpsolver, args);

% lower constraint bounds
lmisol.lbg = lbg;

% retrieve LMI coordinates
sossol = call(obj.lmi2sos,lmisol);

% build polynomial solution
argout = {sossol.x sossol.f sossol.g sossol.lam_x sossol.lam_g};

end
