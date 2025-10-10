classdef (Abstract) ConicSolver < casos.package.solvers.SolverCallback
% Base class for low-level conic (SDP) solvers.
%
% The generic conic problem has the form
%
%   min 1/2*x'*h*x + g'*x
%   st. a*x in Kc, x in Kx, 
%
% with Lagrange multipliers satisfying
%
%   h*x + g + a'*lam_a + lam_x = 0
%   lam_a in -Kc*, lam_x in -Kx*
%
% where Kc and Kx are cones composed of
%   - box constraint [lb ub] (l)
%   - Lorentz (quadratic) cone (q)
%   - rotated Lorentz cone (r)
%   - cone of PSD matrices (s)
%
% and Kc* and Kx* are the dual cones of Kc and Kx, respectively.

properties (GetAccess=protected, SetAccess=private)
    args_in = struct;
    names_out = {'x' 'cost' 'lam_a' 'lam_x'};
end

properties (Access=protected)
    status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
end

properties (Constant, Access=protected)
    conic_options = [casos.package.solvers.SolverCallback.solver_options
        {'hessian_permute', 'Inverse permutation matrix for the Hessian.'}
    ];

    conic_cones = casos.package.Cones([
        casos.package.Cones.LIN
        casos.package.Cones.LOR
        casos.package.Cones.ROT
        casos.package.Cones.PSD
    ]);
end

properties (Abstract, Access=protected)
    fhan;
    ghan;
end

methods (Abstract, Access=protected)
    buildproblem(obj);
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.ConicSolver.conic_options;
    end

    function cones = get_cones
        % Return supported cones.
        cones = casos.package.solvers.ConicSolver.conic_cones;
    end
end

methods
    function obj = ConicSolver(name,conic,varargin)
        obj@casos.package.solvers.SolverCallback(varargin{:});

        % sparsity patterns
        as = conic.a;
        % problem size
        [m,n] = size(as);

        if isfield(conic,'h')
            hs = conic.h;
        else
            hs = casadi.Sparsity(n,n);
        end

        % default options
        if ~isfield(obj.opts,'Kx'), obj.opts.Kx = struct('lin',n); end
        if ~isfield(obj.opts,'Kc'), obj.opts.Kc = struct('lin',m); end
        if ~isfield(obj.opts,'hessian_permute'), obj.opts.hessian_permute = 1; end

        % check cone dimensions
        assert(obj.getnumc(obj.opts.Kx) == n, 'Dimension of Kx must equal to number of variables (%d).', n)
        assert(obj.getnumc(obj.opts.Kc) == m, 'Dimension of Kc must equal to number of constraints (%d).', m)

        % dimensions of linear variables and constraints
        Nl = obj.getdimc(obj.opts.Kx,'lin');
        Ml = obj.getdimc(obj.opts.Kc,'lin');

        % symbolic inputs
        obj.args_in.h      = casadi.MX.sym('h',hs);
        obj.args_in.g      = casadi.MX.sym('g',n);
        obj.args_in.a      = casadi.MX.sym('a',as);
        obj.args_in.lba    = casadi.MX.sym('lba',Ml);
        obj.args_in.uba    = casadi.MX.sym('uba',Ml);
        obj.args_in.cba    = casadi.MX.sym('cba',m-Ml);
        obj.args_in.lbx    = casadi.MX.sym('lbx',Nl);
        obj.args_in.ubx    = casadi.MX.sym('ubx',Nl);
        obj.args_in.cbx    = casadi.MX.sym('cbx',n-Nl);
        obj.args_in.x0     = casadi.MX.sym('x0',n);
        obj.args_in.lam_x0 = casadi.MX.sym('lam_x0',n);
        obj.args_in.lam_a0 = casadi.MX.sym('lam_a0',m);

        % build conic problem
        buildproblem(obj);

        % construct CasADi callback
        construct(obj,name);
    end
end

methods (Access=protected)
    %% Quasi-Static helper functions
    % MOSEK/SCS-style vectorization of the semidefinite cone
    [V,varargout] = sdp_vec(obj,M,varargin);
    [M,varargout] = sdp_mat(obj,V,varargin);

    %% Cone helper functions
    function d = getdimc(obj,K,type)
        % Return cone dimensions of specific type.        
        d = get_dimension(obj.get_cones,K,type);
    end

    function n = getnumc(obj,K,varargin)
        % Return number of variables in cone.
        n = get_length(obj.get_cones,K,varargin{:});
    end
end

end
