function out = str(obj)
% Return string representation for polynomial.

if isempty(obj)
    out = {'[]'};
    return
end

% else
out = cell(size(obj));

% iterate over elements
for i = 1:numel(obj)
    % polynomial terms (monomials)
    terms = cell(1,obj.nterm);

    firstterm = 1;

    for j = 1:obj.nterm
        % coefficient & degree
        cf = obj.coeffs(j,i);
        dg = obj.degmat(j,:);

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
            % zero constant
            continue;
        end

        % remove sign on first term
        if firstterm
            sgn = mdf;
            % reset flag
            firstterm = 0;
        end

        % build monomial
        if sum(dg) == 0
            % constant
            monom = '';
        else
            % monomial
            ms = arrayfun(@exp2str, obj.indets(find(dg)), nonzeros(dg)');
            pd = repmat({'*'},nnz(dg),1);
            % remove leading one
            if is_one(abs(cf))
                scf = '';
                pd{1} = '';
            end

            M = [pd(:)'; ms(:)'];
            monom = [M{:}];
        end

        % combine
        terms{j} = [sgn scf monom];
    end

    % check for zero polynomial
    if firstterm
        out{i} = '0';
    else
        out{i} = [terms{:}];
    end
end

end

function s = exp2str(v,d)
% Convert exponential v^d to string.

    switch d
        case 0, s = {''};
        case 1, s = v;
        otherwise, s = compose('%s^%d',v{:},d);
    end
end