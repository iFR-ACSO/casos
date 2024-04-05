function print_terms(obj)
% Print element-wise monomial terms.

[~,L] = get_degmat(obj);
monom = str_monoms(obj);

% prepare output
out = cell(size(obj));
% preassign chars
out(:) = {'()'};

for i = find(obj)
    % build element-wise monomial terms
    m = join(monom(L(i,:)),',');

    if sum(L(i,:)) > 1
        % embed vector in parentheses
        out(i) = compose('(%s)',m);
    else
        out(i) = m;
    end
end

% print matrix of terms
disp_matrix(obj,'[]',out);

end
