function Z=monomials(vars,deg)
% function Z=monomials(p,deg)
%
% DESCRIPTION
%   Construct list of monomials
%
% INPUTS
%   p: A polynomial, vector of pvars or a non-negative integer.
%   deg: A vector of non-negative integers specifying the degrees
%         of monomials to be included Z.
%
% OUTPUTS
%   Z: lz-by-1 list of monomials
%
%   Z: If p is a polynomial (deg is not specified) then Z will be a
%      lz-by-1 vector of all monomials in p. If p is a vector of pvars
%      then Z will be a lz-by-1 vector of all monomials of the specified
%      degrees in the given pvars. If vars is a non-negative integer then
%      Z will be the lz-by-var degree matrix with each row specifying the
%      degrees of one of the monomials.
%
% SYNTAX
%   Z=monomials(p)
%      If p is a polynomial then Z is the vector of monomials in p.
%   Z=monomials(vars,deg)
%      If vars is a vector of pvars then Z is a vector of all monomials
%      in the variables listed in vars and degrees listed in deg.
%   Z=monomials(nvar,deg)
%      If nvar is a non-negative integer then Z is the degree matrix
%      corresponding to all monomials in nvar variables and degrees deg.
%
% EXAMPLE
%   pvar x1 x2
%   Z1 = monomials([x1;x2],0:2)
%
%   p = x1^2+5*x1*x2-6*x2^3;
%   Z2 = monomials(p)
%
% See also poly2basis

% 4/21/2009: PJS  Initial Coding with single argument call


if nargin==1
    % Single argument -- Construct vector of monomials in p
    p = polynomial(vars);
    if ~isa(p,'polynomial')
        error('p must be a polynomial for Z=monomials(p) syntax');
    end
    
    % Flipping/sorting to return output in lexicographic order
    % (sorted by degree then by alphabetical order)
    degmat = p.degmat;
    tmp = [ sum(degmat,2) fliplr(degmat)];
    degmatsort = sortrows(tmp);
    degmatsort = fliplr(degmatsort(:,2:end));
    
    Nt = size(degmat,1);
    Z = polynomial(speye(Nt),degmatsort,p.varname,[Nt 1]);
    return;
elseif nargin~=2
    error('Syntax: Z=monomials(p,deg)');
end

% Error checking on vars
if isa(vars,'double') && isscalar(vars) && vars>=0 && ...
        floor(vars)==ceil(vars)
    nvar = vars;
    vars = [];
elseif ispvar(vars)
    vars = char(vars);
    nvar = length(vars);
elseif iscellstr(vars)
    vars = vars(:);
    nvar = length(vars);
else
    error('vars must be a scalar positive integer or a vector of pvars');
end

if isa(deg,'double') && ndims(deg)==2 && ...
        all(floor(deg)==ceil(deg)) && all(deg>=0)
    
    % Use recursive call to construct monomials vector
    if nvar==0
        degmat = zeros(1,0);
    elseif nvar==1
        % Base case of recursion
        degmat = deg(:);
    else
        % Recursive call for nvar>1
        deg = deg(:);
        maxd = max(deg);
        
        % Compute all degmats for nvar-1 vars and degs from 0 up to maxd
        degmatcell = cell(maxd+1,1);
        for i1=0:maxd
            degmatcell{i1+1} = monomials(nvar-1,i1);
        end
        
        % Stack on degrees for last variable
        degmat = [];
        for i1 = 1:length(deg)
            degi = deg(i1);
            degmati = [];
            for i3=0:degi
                temp = degmatcell{degi-i3+1};
                temp = [temp repmat(i3,[size(temp,1),1])];
                degmati = [degmati; temp];
            end
            degmat = [degmat;degmati];
        end
    end
    
elseif iscell(deg)
    
    % XXX UNDOCUMENTED
    % New case: deg is a cell array
    % XXX This allows more flexibility in creating a monomials vector,
    % but is this a useful syntax? Or is it easier to just create
    % the list of all monoms with the standard syntax and then prune
    % out ones which are not wanted?
    deg = deg(:);
    ldeg = length(deg);
    if ldeg~=nvar
        error('Length of deg cell array must equal the number of variables');
    end
    
    % Create degree matrix for monomials vector
    degmat = [];
    for i1=1:nvar
        degi = deg{i1};
        degi = degi(:);
        
        if isa(degi,'double') && ndims(degi)==2 && ...
                all(floor(degi)==ceil(degi)) && all(degi>=0)
            if i1==1
                degmat = degi;
            else
                degmat1 = repmat(degmat,[length(degi) 1]);
                degmat2 = repmat(degi',[size(degmat,1) 1]);
                degmat = [degmat1 degmat2(:)];
            end
        else
            error('deg must be a vector or cell array of non-negative integers');
        end
    end
    
else
    error('deg must be a vector or cell array of non-negative integers');
end

% Create output
if isempty(vars)
    % Output degmat if input is number of vars
    Z = degmat;
else
    % Ouput a poly if input is a list of pvars
    if isempty(find(degmat,1))
        Z = polynomial(1);
    else
        ld = size(degmat,1);
        chkval = 0; % skip validity check
        Z = polynomial(speye(ld),degmat,vars,[ld 1],chkval);
    end
end

