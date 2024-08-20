% function opt = sosoptions(Name1,Value1,Name2,Value2,...)
%
% DESCRIPTION
%   Creates an options object for SOSOPT.
%
% INPUTS
%   Name, Value: The options property specified by the character string
%      Name is set to Value. The setable options properties are specified
%      below.  The allowable values for each property are specified in
%      brackets with the default choice in braces.
%
%      -solver: Optimization solver to be used.
%          [{'sedumi'},'sdpam','dsdp','sdpt3','csdp','sdplr','mosek']
%      -form: Formulation for the optimization. [{'image'},'kernel']
%      -simplify: SOS simplification procedure to remove monomials
%         that are not needed in the Gram matrix form. [{'on'}; 'off']
%      -scaling: Scaling of SOS constraints. ['on'; {'off'}]
%      -checkfeas: Check feasibility of solution.  'fast' checks
%         feasibility information returned by the solver. 'full' checks
%         the validity of the Gram matrix decomposition in the output
%         sossol. 'both' does both feasibility checks. 
%         ['off'; {'fast'}; 'full'; 'both']
%      -feastol: Feasibility tolerance used in full feasibility check.
%         [{1e-6}]
%      -solveropts: Structure with options passed directly to the
%         optimization solver.  The default is empty.  The solver
%         display is turned off with this default.
%
% OUTPUT
%   opt: sosoptions object
%
% SYNTAX
%   opt = sosoptions
%     Creates an sosoptions object initialized with default values.
%   opt = sosoptions(Name1,Value1,Name2,Value2,...)
%     Creates an sosoptions object with options specified by Name set
%     to the values specified by Value.
%
% See also sosopt, gsosopt, gsosoptions

% 10/26/2010 PJS  Initial Coding

classdef sosoptions
    
    properties
        solver = 'sedumi';
        form = 'image';
        simplify = 'on';
        scaling = 'off';
        checkfeas = 'fast';
        feastol = 1e-6;
        solveropts = [];
    end
    
    methods
        % Constructor
        function opt = sosoptions(varargin)
            % Check # of input args
            nin = nargin;
            if ceil(nin/2)~=floor(nin/2)
                errstr1 = 'SOSOPTIONS must have an even number of inputs';
                errstr2 = ' with Name/Value pairs specified together.';
                error([errstr1 errstr2]);
            end
            
            % Set Name/Value pairs:
            % Rely on default error if Name is not a public property
            for i1=1:(nin/2)
                Name = varargin{2*(i1-1)+1};
                Value = varargin{2*i1};
                opt.(Name) = Value;
            end
        end
        
        % Set: solver
        % 'setup' is an undocumented call used by gsosopt
        function opt = set.solver(opt,value)
            % XXX Currently working: Sedumi, CSDP, SDPT3, SDPLR, Setup,
            if ischar(value)
                opt.solver = value;
            else
                errstr1 = 'solver can be ''sedumi'', ''sdpam'', ''csdp''';
                errstr2 = ', ''dsdp'', ''sdplr'', ''sdpt3'', or ''mosek''.';
                error([errstr1 errstr2]);
            end
        end
        
        % Set: scaling
        function opt = set.scaling(opt,value)
            AllowableVal = {'on'; 'off'};
            if ischar(value) && any( strcmp(value,AllowableVal) )
                opt.scaling = value;
            else
                error('scaling can be ''on'' or ''off''. ');
            end
        end
        
        % Set: form
        function opt = set.form(opt,value)
            AllowableVal = {'image'; 'kernel'};
            if ischar(value) && any( strcmp(value,AllowableVal) )
                opt.form = value;
            else
                error('form can be ''image'' or ''kernel''. ');
            end
        end
        
        % Set: checkfeas
        function opt = set.checkfeas(opt,value)
            AllowableVal = {'off'; 'fast'; 'full'; 'both'};
            if ischar(value) && any( strcmp(value,AllowableVal) )
                opt.checkfeas = value;
            else
                errstr1 = 'checkfeas can be ''off'', ''fast'',';
                errstr2 = ' ''full'', or ''both''.';
                error([errstr1 errstr2]);
            end            
        end
        
        % Set: feastol
        function opt = set.feastol(opt,value)
            if isscalar(value) && isa(value,'double') && value>0
                opt.feastol = value;
            else
                error('feastol must be a positive, scalar number.');
            end
        end
        
        % Set: simplify
        function opt = set.simplify(opt,value)
            AllowableVal = {'on'; 'off'};
            if ischar(value) && any( strcmp(value,AllowableVal) )
                opt.simplify = value;
            else
                error('simplify can be ''on'' or ''off''. ');
            end
        end
        
        % Set: solveropts
        function opt = set.solveropts(opt,value)
            opt.solveropts = value;
        end
        
    end % methods
end % classdef