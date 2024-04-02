function B = cleanpoly(p,tol,deg)
% Remove coefficients from nonsymbolic polynomials.

% check inputs
if nargin<3
    deg = [];
end

assert(~is_symexpr(p), 'Coefficients must not be symbolic expressions.')

% copy coefficients
coeffs = casadi.SX(p.coeffs);

% sparse zero
sp_zero = casadi.SX(1,1);

if ~isempty(tol)
    % remove coefficients below tolerance
    assert(isdouble(tol) && isscalar(tol) && tol >= 0, 'Tolerance must be a positive scalar or is empty.')

    % polynomial coefficients are casadi SX; we need numerical values
    idx = find(abs(full(casadi.DM(coeffs))) < tol);

    if ~isempty(idx)
        % remove coefficients that are below tolerance
        coeffs(idx) = sp_zero;
    end
end

if ~isempty(deg)
    % find all parts that shall be retained with certain degrees
    if isa(deg,'double') &&  all(floor(deg)==ceil(deg)) && all(deg>=0)
        % degree of each term
        pdeg = sum(p.degmat,2);

        % find all terms we don't want to retain
        idx = find(~ismember(pdeg,deg));
        
        % remove the coefficients
        coeffs(idx,:) = sp_zero; %#ok<FNDSB> 

    % find all parts that shall be retained with certain degrees of a variable
    elseif iscell(deg)  && size(deg,2)== 2
        warning('This use of cleanpoly is deprecated.')

        % deg is a cell array of pvars and vectors of non-neg integers
        for i1=1:size(deg,1)
            % check if variable is part of indeterminates
            var = deg{i1,1};

            vidx = find(strcmp(p.indets,var));
            
            if ~isempty(vidx)
                
                degi = deg{i1,2};
                
                % check which degrees shall be retained
                assert(isa(degi,'double') && size(degi,2) == 2 && all(floor(degi) == ceil(degi)) && all(degi>=0), ...
                    'Second column of cell array deg must contain a vector of non-negative integers.' ...
                );
                
                % get entries in degmat of the desired indeterminates we want to retain
                pdeg = p.degmat(:,vidx);

                idx = arrayfun(@(x) find(pdeg == x), unique(degi), 'UniformOutput', false);
                idx = vertcat(idx{:});

                % get entries we don't want to retain
                idx = setdiff(1:length(pdeg),idx);
                coeffs(idx,:) = sp_zero;
            end
        end

    else
        error('Third input must be a vector of non-negative integers.');
    end
end

% new polynomial
B = casos.PS;

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,p.degmat,p.indets);

B.coeffs = coeffs;
B.degmat = degmat; 
B.indets = indets;
B.matdim = p.matdim;

end
