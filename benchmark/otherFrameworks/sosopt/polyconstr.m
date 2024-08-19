% function pco = polyconstr(LeftSide,RightSide,RelOp);
%
% DESCRIPTION
%   Creates a polynomial constraint object for SOSOPT and GSOSOPT.
%   The relational operators >=, <=, and == are overloaded so that p>=q,
%   p<=q, and p==q will generate a polynomial constraint object from
%   polynomial vectors p and q.
%
% INPUTS
%   LeftSide: N-by-1 polynomial vector on left side of relational operator
%   RightSide: N-by-1 polynomial vector on right side of relational operator
%   RelOp: Relational operator. RelOp can be <=, >=, or ==.
%
% OUTPUT
%   pco: polynomial constraint object
%
% SYNTAX
%   pco = polyconstr
%     Returns an empty polyconstr object
%   pco = polyconstr(LeftSide,RightSide,RelOp)
%     Creates a polynomial constraint object.
%   pco = p>=q
%   pco = p<=q
%   pco = p==q
%     Creates a polynomial constraint object. If p is a scalar and q is a
%     vector then p will be expanded to the same dimension as q. Similarly,
%     if q is a scalar and p is a vector then q will be expanded. The
%     data for the polynomial constraint can be accessed by pco.LeftSide,
%     pco.RightSide, and pco.RelOp.  In addition, the constraint is
%     converted to an equivalent one-sided form, pco.OneSide, where the
%     relation is either == or >= and the right side is zero. The
%     constraint is displayed in its one-sided form.
%
% See also sosopt, gsosopt

% 10/22/2010  PJS  Initial Coding (N-by-1 array with Mi-by-1 polys)
% 11/28/2010  PJS  Switched object to 1-by-1 array with N-by-1 polys

classdef polyconstr
    
    % XXX PJS: Is there any reason to keep LeftSide and RightSide? I'm
    % wondering if it would be simpler to have polyconstr inherit all
    % polynomial properties and simply add one field 'RelOp' with 
    % dimensions equal to the polynomial dimensions. Maybe this
    % would cause confusion because only a subset of polynomial
    % methods would work on polyconstr. Moreover, the polyconstr methods
    % would have to do something different than the corresponding
    % polynomial methods.
    
    properties
        LeftSide = polynomial([]);      % polynomial
        RightSide = polynomial([]);     % polynomial
        RelOp = cell(0);                % >=, <=, ==
    end
    properties (Dependent = true, SetAccess = private)
        OneSide                         % Converted to s(x)==0 or s(x)>=0
    end

    methods
        % Constructor
        function pco = polyconstr(LeftSide,RightSide,RelOp)
            if nargin==3
                try
                    LeftSide =  polynomial(LeftSide(:));
                    RightSide = polynomial(RightSide(:));
                catch
                    error('Both sides of inequality must be polynomials.');
                end
                if ~isempty(RelOp) && ~any( strcmp(RelOp,{'<=','>=','=='}) )
                    error('RelOp must be <=, >= or ==');
                end
                
                % Check vector dimensions and perform scalar expansion
                nl = length(LeftSide);
                nr = length(RightSide);
                if nl>1 && nr==1
                    RightSide = repmat(RightSide,[nl 1]);
                elseif nl==1 && nr>1
                    LeftSide = repmat(LeftSide,[nr 1]);
                    nl = nr;
                elseif nl~=nr
                    error('Vector dimensions of polynomials must agree');
                end
                
                % Fill in object fields
                pco.LeftSide = LeftSide;
                pco.RightSide = RightSide;
                pco.RelOp = repmat({RelOp},[nl 1]);
            elseif nargin~=0
                errstr1 = 'Invalid syntax for the "polyconstr" command.';
                errstr2 = ' Type "help polyconstr" for more information.';
                error([errstr1 errstr2]);
            end
        end % polyconstr constructor
        
        % get.OneSide
        function s = get.OneSide(p)
            % Convert to one-sided constraint: s(x)==0 or s(x)>=0
            s = p.LeftSide - p.RightSide;
            idx = find( strcmp('<=',p.RelOp) );
            s(idx) = -s(idx);
        end % get.OneSide
        
        % XXX PJS--Include get functions for polynomial fields and
        % have them simply get the corresponding fields in OneSide?
        % This would increase the connection between polyconstr and
        % polynomial objects.
        %         % get.varname
        %         function out = get.varname(p)
        %             out = unique( [p.LeftSide.varname; p.RightSide.varname] );
        %         end % get.OneSide
        %         % get.coefficient
        %         function out = get.coefficient(p)
        %             out = get.OneSide(p);
        %             out = out.coefficient
        %         end % get.OneSide
        %         % get.degmat
        %         % get.matdim
        %         % nterms;
        %         % nvars;
        %         % maxdeg;
        %         % mindeg;        
        
        % vertcat
        function out = vertcat(varargin)
            if nargin==1
                out = varargin{1};
            else
                v1 = varargin{1};
                if isempty(v1)
                    v1 = polyconstr;
                end
                v2 = varargin{2};
                if isempty(v2)
                    v2 = polyconstr;
                end
                out = polyconstr;
                out.LeftSide = vertcat(v1.LeftSide,v2.LeftSide);
                out.RightSide = vertcat(v1.RightSide,v2.RightSide);
                out.RelOp = vertcat(v1.RelOp,v2.RelOp);
            end
            if nargin>2
                out = vertcat(out,varargin{3:end});
            end
            %             if all(size(out)>1)
            %                 error(['Constraints may be stacked horizontally or '...
            %                     'vertically but not both.']);
            %             end
        end % vertcat
        
        % horzcat
        function out = horzcat(varargin)
            error('Constraints must be stacked vertically.');
            
            %             if nargin==1
            %                 out = varargin{1};
            %             else
            %                 v1 = varargin{1};
            %                 v2 = varargin{2};
            %                 out = polyconstr;
            %                 out.LeftSide = horzcat(v1.LeftSide,v2.LeftSide);
            %                 out.RightSide = horzcat(v1.RightSide,v2.RightSide);
            %                 out.RelOp = horzcat(v1.RelOp,v2.RelOp);
            %             end
            %             if nargin>2
            %                 out = horzcat(out,varargin{3:end});
            %             end
            %             if all(size(out)>1)
            %                 error(['Constraints must be stacked horizontally or '...
            %                     'vertically but not both.']);
            %             end
        end % horzcat
        
        % size
        function varargout = size(a,dim)
            if nargin ==1
                out = size(a.LeftSide);
                if nargout==0 || nargout==1
                    varargout{1}=out;
                elseif nargout==2
                    varargout{1} = out(1);
                    varargout{2} = out(2);
                end
            elseif nargin == 2
                out = size(a.LeftSide,dim);
                varargout{1} = out;
            end
        end % size
        
        % length
        function b = length(a)
            b = length(a.LeftSide);
        end
        
        % isempty
        function b = isempty(a)
            b = isempty(a.LeftSide);
        end
        
        % end
        function out = end(a,k,n)
            if n ==1
                sza = size(a);
                out = sza(1)*sza(2);
            elseif n==2
                out = size(a,k);
            end
        end
        
        % subsref
        function b = subsref(a,L)
            b = a;
            switch L(1).type
                case '.'
                    b = a.(L(1).subs);
                case '()'
                    b.LeftSide = subsref(b.LeftSide,L(1));
                    b.RightSide = subsref(b.RightSide,L(1));
                    b.RelOp = subsref(b.RelOp,L(1));
                case '{}'
                    error(['{}- like subsreference is not supported '...
                        'for polyconstr objects.']);
            end
            if length(L)>1
                b = subsref(b,L(2:end));
            end
        end
        
        % subsasgn
        function b = subsasgn(a,L,RHS)
            if isempty(a)
                a = polyconstr;
            end
            switch L(1).type
                case '.'
                    if length(L) == 1
                        temp = RHS;
                    else
                        temp = subsref(a,L(1));
                        temp = subsasgn(temp,L(2:end),RHS);
                    end
                    b = a;
                    b.(L(1).subs) = temp;
                case '()'
                    % Peform all subsasgn but L(1)
                    if length(L)==1
                        temp = RHS;
                    else
                        temp = subsref(a,L(1));
                        temp = subsasgn(temp,L(2:end),RHS);
                    end
                    
                    b = a;
                    b.LeftSide = subsasgn(b.LeftSide,L(1),temp.LeftSide);
                    b.RightSide = subsasgn(b.RightSide,L(1),temp.RightSide);
                    b.RelOp = subsasgn(b.RelOp,L(1),temp.RelOp);
                    
                    b.LeftSide = b.LeftSide(:);
                    b.RightSide = b.RightSide(:);
                    b.RelOp = b.RelOp(:);
                case '{}'
                    error(['{}- like subsassign is not supported '...
                        'for polyconstr objects.']);
            end
        end
        
        % display
        function display(p)
            cflag = isequal(get(0,'FormatSpacing'),'compact');
            Np = prod(size(p.LeftSide));
            if ~cflag
                disp(' ');
            end
            if Np == 0
                disp([inputname(1) ' = ']);
                display(' Empty polyconstr object.');
            elseif Np==1
                RelOp = p.RelOp;
                
                eval([inputname(1) '= p.OneSide;']);
                eval(['display(' inputname(1) ');']);
                if strcmp(RelOp,'==')
                    display('  == 0');
                else
                    display('  >= 0');
                end
            else
                disp([inputname(1) ' = ']);
                display(['  polyconstr object with ' int2str(Np) ...
                    ' constraints.']);
            end
            if ~cflag
                disp(' ');
            end
        end % display
        
        % subs
        function B = subs(p,varargin)
            B = p;
            B.LeftSide = subs(B.LeftSide,varargin{:});
            B.RightSide = subs(B.RightSide,varargin{:});
        end % subs
        
    end % methods
end % classdef