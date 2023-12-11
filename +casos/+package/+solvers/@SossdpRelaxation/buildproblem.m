function buildproblem(obj,solver,sos)
% Build SDP problem from SOS relaxation.

opts = obj.opts;

% problem size
n = length(sos.x);
m = length(sos.g);

% get cone dimensions
if isfield(opts.Kx,'l'), Nl = opts.Kx.l; else, Nl = 0; end
if isfield(opts.Kx,'s'), Ns = opts.Kx.s; else, Ns = 0; end
if isfield(opts.Kc,'l'), Ml = opts.Kc.l; else, Ml = 0; end
if isfield(opts.Kc,'s'), Ms = opts.Kc.s; else, Ms = 0; end

assert(n == (Nl + Ns), 'Dimension of Kx must be equal to number of variables (%d).', n);
assert(m == (Ml + Ms), 'Dimension of Kc must be equal to number of constraints (%d).', m)

% select sum-of-squares variables and constraints
Is = [false(Nl,1); true(Ns,1)];
Js = [false(Ml,1); true(Ms,1)];

% obtain Gram basis and matrix decision variables
[Qgram_x,Zgram_x,Ksdp_x_s] = grammatrix(sos.x,Is);
% obtain Gram basis for sum-of-squares constraints
[Zgram_g,Ksdp_g_s] = grambasis(sos.g,Js);

assert(length(Qgram_x) == size(Zgram_x,1), 'Unable to find Gram matrix of decision variables.')

% matrix decision variables for constraints
Qgram_g = casadi.SX.sym('Q',sum(Ksdp_g_s.^2),1);

% replace sum-of-squares constraints by equality
gdiff = (sos.g - [zeros(Ml,1); casos.PS(Zgram_g,Qgram_g)]);
% handle (new) equality constraints
[Qdiff_g,Zdiff_g] = poly2basis(gdiff);

% handle linear decision variables
[Qlin_x,Zlin_x] = poly2basis(sos.x,~Is);
% handle parameter
[Qlin_p,Zlin_p] = poly2basis(sos.p);

% get cost function
[Qlin_f,Zlin_f] = poly2basis(sos.f);

% build SDP problem
sdp.x = [Qlin_x; Qgram_x; Qgram_g];
sdp.f = Qlin_f;
sdp.g = Qdiff_g;
sdp.p = Qlin_p;
% SDP options
sdpopt = opts.sdpsol_options;
sdpopt.Kx = struct('l', numel(Qlin_x), 's', [Ksdp_x_s; Ksdp_g_s]);
sdpopt.Kc = struct('l', numel(Qdiff_g));

% initialize SDP solver
obj.sdpsolver = casos.sdpsol('SDP',solver,sdp,sdpopt);

% store basis
obj.monom_xl = Zlin_x;
obj.monom_xs = basis(sos.x,Is);
obj.monom_p  = Zlin_p;
obj.monom_f  = Zlin_f;
obj.monom_gl = basis(gdiff,~Js);
obj.monom_gs = basis(gdiff, Js);
% gram basis
obj.gram_x = Zgram_x;
obj.gram_g = Zgram_g;
% output basis
obj.basis_x_out = blkdiag(obj.monom_xl, obj.gram_x);
obj.basis_g_out = blkdiag(obj.monom_gl, obj.gram_g);

end
