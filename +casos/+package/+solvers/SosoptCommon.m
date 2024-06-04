classdef (Abstract) SosoptCommon < casos.package.solvers.SolverCommon & casos.package.functions.FunctionInterface
% Common superclass for sum-of-squares optimization problems.
%
% The generic sum-of-squares problem has the form
%
%   min f(x,p) s.t. g(x,p) in Kc, x in Kx
%
% with Lagrange multipliers satisfying
%
%   df(x,p) + dg(x,p)[lam_g] + lam_x = 0
%   lam_g in -Kc*, lam_x in -Kx*
%
% where Kc and Kx are cones composed of
%   - coefficient-wise inequalities (l)
%   - sum-of-squares polynomials (s)
%
% and Kc* and Kx* are the dual cones of Kc and Kx, respectively.
%

properties (Constant,Access=protected)
    sosopt_options = casos.package.solvers.SolverCommon.solver_options;

    sosopt_cones = casos.package.Cones([
        casos.package.Cones.LIN
        casos.package.Cones.SOS
    ]);

    name_i = {'x0' 'p' 'lbx' 'ubx' 'cbx' 'lbg' 'ubg' 'cbg' 'lam_x0' 'lam_g0'};
    name_o = {'x' 'f' 'g' 'lam_x' 'lam_g'};
end

properties (Access=protected)
    sparsity_xl;
    sparsity_xs;
    sparsity_p;
    sparsity_f;
    sparsity_gl;
    sparsity_gs;
end

properties (Access=protected,Dependent)
    sparsity_x;
    sparsity_g;
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SosoptCommon.sosopt_options;
    end

    function cones = get_cones
        % Return supported cones.
        cones = casos.package.solvers.SosoptCommon.sosopt_cones;
    end
end

methods
    function [obj,sos] = SosoptCommon(name,sos,varargin)
        obj@casos.package.functions.FunctionInterface(name);
        obj@casos.package.solvers.SolverCommon(varargin{:});

        % problem size
        if isfield(sos,'x')
            n = length(sos.x);
        else
            n = 0;
        end
        if isfield(sos,'g')
            m = length(sos.g);
        else
            m = 0;
        end

        % default options
        if ~isfield(obj.opts,'Kx'), obj.opts.Kx = struct('lin',n); end
        if ~isfield(obj.opts,'Kc'), obj.opts.Kc = struct('lin',m); end
    end

    %% Getter
    function z = get.sparsity_x(obj)
        % Sparsity of decision variables.
        z = vertcat(obj.sparsity_xl,obj.sparsity_xs);
    end

    function z = get.sparsity_g(obj)
        % Sparsity of constraints.
        z = vertcat(obj.sparsity_gl,obj.sparsity_gs);
    end

    function n = get_n_in(obj)
        % Number of inputs.
        n = length(obj.name_i);
    end

    function str = get_name_in(obj,i)
        % Name of inputs.
        if nargin > 1
            str = obj.name_i{i+1};
        end
    end

    function z = get_sparsity_in(obj,idx)
        % Sparsity of inputs.
        switch (idx)
            case 0, z = obj.sparsity_x;
            case 1, z = obj.sparsity_p;
            case {2 3}, z = obj.sparsity_xl;
            case 4, z = obj.sparsity_xs;
            case {5 6}, z = obj.sparsity_gl;
            case 7, z = obj.sparsity_gs;
            case 8, z = obj.sparsity_x;
            case 9, z = obj.sparsity_g;
            otherwise, error('Index out of bound (%d).',idx);
        end
    end

    function val = get_default_in(~,i)
        % Defaults of inputs.
        switch (i)
            case {2 5}, val = -inf;
            case {3 6}, val = +inf;
            otherwise, val = 0;
        end
    end

    function idx = get_index_in(obj,str)
        % Index of inputs.
        [tf,ii] = ismember(str,obj.name_i);

        assert(tf, 'Could not find entry "%s".', str);

        % zero-based index
        idx = ii - 1;
    end

    function n = get_n_out(obj)
        % Number of outputs.
        n = length(obj.name_o);
    end

    function str = get_name_out(obj,idx)
        % Name of outputs.
        if nargin > 1
            str = obj.name_o{idx+1};
        end
    end

    function z = get_sparsity_out(obj,idx)
        % Sparsity of outputs.
        switch (idx)
            case 0, z = obj.sparsity_x;
            case 1, z = obj.sparsity_f;
            case 2, z = obj.sparsity_g;
            case 3, z = obj.sparsity_x;
            case 4, z = obj.sparsity_g;
            otherwise, error('Index out of bound (%d).',idx);
        end
    end

    function idx = get_index_out(obj,str)
        % Index of outputs.
        [tf,ii] = ismember(str,obj.name_o);

        assert(tf, 'Could not find entry "%s".', str);

        % zero-based index
        idx = ii - 1;
    end
end

end
