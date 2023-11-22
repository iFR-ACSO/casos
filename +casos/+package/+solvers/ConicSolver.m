classdef (Abstract) ConicSolver < casos.package.solvers.SolverCallback
% Base class for low-level conic (SDP) solvers.
%
% The generic conic problem has the form
%
%   min 1/2*x'*h*x + g'*x
%   st. a*x in Ka, x in Kx, 
%
% where Ka and Kx are cones composed of
%   - box constraint [lb ub] (l)
%   - Lorentz (quadratic) cone (q)
%   - rotated Lorentz cone (r)
%   - cone of PSD matrices (s)
%

properties (GetAccess=protected, SetAccess=private)
    args_in = struct;
    names_out = {'x' 'cost' 'lam_a' 'lam_x'};
end

properties (Access=protected)
    status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
end

properties (Constant, Access=protected)
    conic_options = [casos.package.solvers.SolverCallback.solver_options
        {'Kx', 'Cone description for state constraints.'
         'Ka', 'Cone description for constraint function.'}
    ];
end

properties (Abstract, Access=protected)
    fhan;
    ghan;
end

methods (Abstract, Access=protected)
    buildproblem(obj);
end

methods (Static, Access=protected)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.ConicSolver.conic_options;
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
        if ~isfield(obj.opts,'Kx'), obj.opts.Kx = struct('l',n); end
        if ~isfield(obj.opts,'Ka'), obj.opts.Ka = struct('l',m); end

        % check cone dimensions
        assert(sum(cellfun(@(fn) obj.getnumc(obj.opts.Kx,fn), fieldnames(obj.opts.Kx))) == n, 'Dimension of Kx must equal to number of variables (%d).', n)
        assert(sum(cellfun(@(fn) obj.getnumc(obj.opts.Ka,fn), fieldnames(obj.opts.Ka))) == m, 'Dimension of Ka must equal to number of constraints (%d).', m)

        % dimensions of linear variables and constraints
        Nl = casos.package.solvers.ConicSolver.getdimc(obj.opts.Kx,'l');
        Ml = casos.package.solvers.ConicSolver.getdimc(obj.opts.Ka,'l');

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

methods (Static, Access=protected)
    %% Static helper functions
    function N = getdimc(K,type)
        % Return cone dimensions of specific type.        
        if isfield(K,type)
            N = K.(type);
        elseif strcmp(type,'l')
            N = 0;
        else
            N = [];
        end
    end

    function n = getnumc(K,type)
        % Return number of variables in cone.
        if ~isfield(K,type)
            n = 0;
            return
        end

        % else
        switch (type)
            case 'l', n = K.(type);
            case 's', n = sum(K.(type).^2);
            otherwise, n = sum(K.(type));
        end
    end
end

end
