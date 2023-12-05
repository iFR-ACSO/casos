classdef (Abstract) SosoptCommon < casos.package.functions.FunctionInterface
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
    sosopt_options = [casos.package.functions.FunctionInterface.options
        {'Kx', 'Cone description for state constraints.'
         'Kc', 'Cone description for constraint function.'}
    ];

    name_i = {'x0' 'p' 'lbx' 'ubx' 'cbx' 'lbg' 'ubg' 'cbg' 'lam_x0' 'lam_g0'};
    name_o = {'x' 'f' 'g' 'lam_x' 'lam_g'};
end

properties (Access=protected)
    monom_xl;
    monom_xs;
    monom_p;
    monom_f;
    monom_gl;
    monom_gs;
end

properties (Access=protected,Dependent)
    monom_x;
    monom_g;
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SosoptCommon.sosopt_options;
    end
end

methods
    function [obj,sos] = SosoptCommon(name,sos,varargin)
        obj@casos.package.functions.FunctionInterface(name,varargin{:});

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
        if ~isfield(obj.opts,'Kx'), obj.opts.Kx = struct('l',n); end
        if ~isfield(obj.opts,'Kc'), obj.opts.Kc = struct('l',m); end
        if ~isfield(obj.opts,'sdpsol_options'), obj.opts.sdpsol_options = struct; end
    end

    %% Getter
    function z = get.monom_x(obj)
        % Monomials of decision variables.
        z = blkdiag(obj.monom_xl,obj.monom_xs);
    end

    function z = get.monom_g(obj)
        % Monomials of constraints.
        z = blkdiag(obj.monom_gl,obj.monom_gs);
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

    function z = get_monomials_in(obj,idx)
        % Monomials of inputs.
        switch (idx)
            case 0, z = obj.monom_x;
            case 1, z = obj.monom_p;
            case {2 3}, z = obj.monom_xl;
            case 4, z = obj.monom_xs;
            case {5 6}, z = obj.monom_gl;
            case 7, z = obj.monom_gs;
            case 8, z = obj.monom_x;
            case 9, z = obj.monom_g;
            otherwise, error('Index out of bound (%d).',idx);
        end
    end

    function sz = get_size_in(obj,i)
        % Size of inputs.
        sz = [size(get_monomials_in(obj,i),2) 1];
    end

    function val = get_default_in(~,i)
        % Defaults of inputs.
        switch (i)
            case {2 5}, val = casos.PS(-inf);
            case {3 6}, val = casos.PS(+inf);
            otherwise, val = casos.PS(0);
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

    function z = get_monomials_out(obj,idx)
        % Monomials of outputs.
        switch (idx)
            case 0, z = obj.monom_x;
            case 1, z = obj.monom_f;
            case 2, z = obj.monom_g;
            case 3, z = obj.monom_x;
            case 4, z = obj.monom_g;
            otherwise, error('Index out of bound (%d).',idx);
        end
    end

    function sz = get_size_out(obj,i)
        % Size of outputs.
        sz = [size(get_monomials_out(obj,i),2) 1];
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
