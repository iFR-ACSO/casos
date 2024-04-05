function out = str_monoms(obj)
% Return string representation for monomials.

out = cell(1,obj.nterm);

for it = 1:obj.nterm
    % degrees
    dg = obj.degmat(it,:);

    if sum(dg) == 0
        % constant
        out{it} = '1';
    else
        % monomial
        ms = arrayfun(@exp2str, obj.indets(find(dg)), nonzeros(dg)');

        out(it) = join(ms,'*');
    end
end

end

function s = exp2str(v,d)
% Convert exponential v^d to string.

    v = str(v);
    switch d
        case 0, s = {''};
        case 1, s = v;
        otherwise, s = compose('%s^%d',v{:},d);
    end
end
