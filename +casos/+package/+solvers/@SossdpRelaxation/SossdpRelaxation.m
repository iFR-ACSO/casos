classdef SossdpRelaxation < casos.package.functions.FunctionInterface
% Solve sum-of-squares problem by relaxation to SDP.

properties (Access=private,Constant)
    name_i = {'x0' 'p' 'lbx' 'ubx' 'cbx' 'lbg' 'ubg' 'cbg' 'lam_x0' 'lam_g0'};
    name_o = {'x' 'f' 'g' 'lam_x' 'lam_g'};
end

properties (Access=private)
    sdpsolver;

    sosopt;

    monom_xl;
    monom_xs;
    monom_p;
    monom_f;
    monom_gl;
    monom_gs;

    gram_x;
    gram_g;
end

properties (Access=private,Dependent)
    monom_x;
    monom_g;
end

properties (SetAccess=protected)
    class_name = 'SossdpRelaxation';
end

methods
    out = call(obj,in);

    function obj = SossdpRelaxation(name,solver,sos,varargin)
        obj@casos.package.functions.FunctionInterface(name);

        % parameter
        if ~isfield(sos,'p')
            sos.p = casos.PS;
        end
        % objective
        if ~isfield(sos,'f')
            % feasibility problem
            sos.f = casos.PS(0);
        end

        % problem size
        n = length(sos.x);
        m = length(sos.g);

        % default options
        obj.sosopt.Kx = struct('l',n);
        obj.sosopt.Kg = struct('s',m);
        obj.sosopt.sdpoptions = struct;
        obj.sosopt.error_on_fail = true;

        % parse options
        [obj.sosopt,~] = casos.package.solvers.parse_options(obj.sosopt,varargin{:});

        % build SDP problem
        obj = buildproblem(obj,solver,sos);
    end
end

methods
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

    function s = get_stats(obj)
        % Return stats.
        s = obj.sdpsolver.stats;
    end
end

end
