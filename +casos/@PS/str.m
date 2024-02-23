function out = str(obj)
% Return string representation for polynomial.

if isempty(obj)
    out = {'[]'};
    return
end

% else
out = cell(size(obj));

% coefficient sparsity
S = sparsity(obj.coeffs);

% monomials
monom = cell(1,obj.nterm);

% polynomial terms
terms = cell(size(S));
% preassign chars
terms(:) = {''};

% check whether element has been visited yet
firstterm = true(1,obj.numel);

% iterate over nonzero coefficients
for ic = find(S)
    % get term and element indices
    [it,ie] = ind2sub(size(S),ic);

    % coefficient & degree
    cf = obj.coeffs(ic);
    dg = obj.degmat(it,:);

    % string representation of coefficient and sign
    mdf = '';
    if is_symbolic(cf)
        % symbolic coefficient
        scf = sprintf('(%s)', str(cf));
        sgn = ' + ';
    elseif ~is_constant(cf)
        % symbolic expression
        scf = str(cf);
        sgn = ' + ';
    elseif is_one(cf > 0)
        % positive (constant) coefficient
        scf = str(cf);
        sgn = ' + ';
    elseif is_one(cf < 0)
        %  negative (constant) coefficient
        scf = str(abs(cf));
        sgn = ' - ';
        mdf = '-';
    else
        % zero constant (nonsparse zero)
        scf = '0';
        sgn = ' + ';
    end

    % remove sign on first term
    if firstterm(ie)
        sgn = mdf;
        % mark element as visited
        firstterm(ie) = 0;
    end

    % build monomial
    if ~isempty(monom{it})
        % nothing to do
    elseif sum(dg) == 0
        % constant
        monom{it} = '';
    else
        % monomial
        ms = arrayfun(@exp2str, obj.indets(find(dg)), nonzeros(dg)');

        monom{it} = strjoin(ms,'*');
    end

    % remove leading one
    if sum(dg) == 0
        % constant
        pd = '';
    elseif is_one(abs(cf))
        scf = '';
        pd = '';
    else
        pd = '*';
    end

    % combine
    terms{ic} = [sgn scf pd monom{it}];
end

% combine terms 
out(~firstterm) = join(terms(:,~firstterm)', '', 2);
% assign sparse zeros
out(firstterm) = {'00'};

end

function s = exp2str(v,d)
% Convert exponential v^d to string.

    switch d
        case 0, s = {''};
        case 1, s = v;
        otherwise, s = compose('%s^%d',v{:},d);
    end
end