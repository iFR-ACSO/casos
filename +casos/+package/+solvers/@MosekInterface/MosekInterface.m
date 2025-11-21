classdef (Sealed) MosekInterface < casos.package.solvers.ConicSolver
    % Interface for conic solver MOSEK.

    properties (Access=protected)
        fhan;
        ghan;
        barv;
        cone;
    end

    properties (Access=private)
        info = struct;
    end

    properties (Constant, Access=protected)
        mosek_options = [casos.package.solvers.ConicSolver.conic_options
            {'mosek_param', 'Parameters to be passed to MOSEK.'
            'cholesky_method','Parameter that defines how cholseky is computed (symbollically or numerically online) '
            'mosek_echo',  'Verbosity level passed to MOSEK (default: 0).'
            'augmented_check', 'Struct with lowered tolerances to find an acceptable solution.'}
            ];
    end

    methods (Static)
        function options = get_options
            % Return static options.
            options = casos.package.solvers.MosekInterface.mosek_options;
        end
    end

    methods
        %     init(obj);
        argout = eval(obj,argin);

        function obj = MosekInterface(name,conic,varargin)
            % Construct MOSEK interface.
            obj@casos.package.solvers.ConicSolver(name,conic,varargin{:});


            % default options
            if ~isfield(obj.opts,'mosek_param'), obj.opts.mosek_param = struct; end
            if ~isfield(obj.opts,'mosek_echo'), obj.opts.mosek_echo = 0; end
            if ~isfield(obj.opts,'augmented_check')
                % default values for augmented solution check, in case mosek
                % returns unknown as solution status
                obj.opts.augmented_check.optMeas     = 0.5;  % Compute optimality measure from Mosek's solver status
                % It should converge to +1 to be optimal% Ensure relative gap is less than x percen
                obj.opts.augmented_check.feasTol     = 1e-6; % relaxed threshold
                obj.opts.augmented_check.relativeGap = 0.05; % Ensure relative gap is less than X percent (default 5%)
                obj.opts.augmented_check.maxNorm     = 1e10; % Avoid using solutions with extremely large norms (ill-conditioning)



            end

        end

        function s = stats(obj)
            % Return stats.
            s = obj.info;
            s = addfield(obj.status,s);
        end

        function sp = get_sparsity_in(obj,i)
            % Return sparsity pattern.
            if (i == 0)
                % Hessian pattern
                sp = sparsity(obj.args_in.h);

            else
                sp = get_sparsity_in@casos.package.solvers.ConicSolver(obj,i);
            end
        end
    end

    methods (Access=protected)
        buildproblem(obj);
    end

end
