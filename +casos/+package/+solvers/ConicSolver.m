classdef (Abstract) ConicSolver < casadi.Callback
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
    sdpopt;

    args_in = struct;
    names_out = {'x' 'cost' 'lam_a' 'lam_x'};
end

properties (Abstract, Access=protected)
    fhan;
    ghan;
end

methods (Abstract, Access=protected)
    buildproblem(obj);
end

methods
    function obj = ConicSolver(name,conic,opts)
        obj@casadi.Callback;

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
        obj.sdpopt.Kx = struct('l',n);
        obj.sdpopt.Ka = struct('l',m);
        obj.sdpopt.solveroptions = [];
        obj.sdpopt.error_on_fail = true;

        % parse options
        if nargin < 3
            opts = struct;
        end
        for fn=intersect(fieldnames(obj.sdpopt),fieldnames(opts))
            if isempty(fn), continue; end
            obj.sdpopt.(fn{:}) = opts.(fn{:});
            opts = rmfield(opts,fn);
        end

        % check cone dimensions
        assert(sum(cellfun(@(fn) obj.getnumc(obj.sdpopt.Kx,fn), fieldnames(obj.sdpopt.Kx))) == n, 'Dimension of Kx must equal to number of variables (%d).', n)
        assert(sum(cellfun(@(fn) obj.getnumc(obj.sdpopt.Ka,fn), fieldnames(obj.sdpopt.Ka))) == m, 'Dimension of Ka must equal to number of constraints (%d).', m)

        % dimensions of linear variables and constraints
        Nl = casos.package.solvers.ConicSolver.getdimc(obj.sdpopt.Kx,'l');
        Ml = casos.package.solvers.ConicSolver.getdimc(obj.sdpopt.Ka,'l');

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
        construct(obj,name,opts);
    end

    %% Common Callback interface
    function n = get_n_in(obj)
        % Return number of inputs.
        n = n_in(obj.fhan);
    end

    function n = get_n_out(obj)
        % Return number of outputs.
        n = n_out(obj.ghan);
    end

    function str = get_name_in(obj,i)
        % Return names of input arguments.
        str = name_in(obj.fhan,i);
    end

    function str = get_name_out(obj,i)
        % Return names of output arguments.
        str = name_out(obj.ghan,i);
    end

    function sp = get_sparsity_in(obj,i)
        % Return sparsity of input arguments.
        sp = sparsity_in(obj.fhan,i);
    end

    function sp = get_sparsity_out(obj,i)
        % Return sparsity of output arguments.
        sp = sparsity_out(obj.ghan,i);
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