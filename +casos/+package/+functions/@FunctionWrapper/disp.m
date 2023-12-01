function disp(obj)
% Display function signature.

args_i = cell(1,obj.n_in);
args_o = cell(1,obj.n_out);

for idx=0:obj.n_in-1
    args_i{idx+1} = get_signature(name_in(obj,idx),monomials_in(obj,idx),size_in(obj,idx));
end
for idx=0:obj.n_out-1
    args_o{idx+1} = get_signature(name_out(obj,idx),monomials_out(obj,idx),size_out(obj,idx));
end

% separate arguments by comma
args_i(2,1:end-1) = {','};
args_o(2,1:end-1) = {','};

fprintf('%s:(%s)->(%s) %s\n', obj.name, [args_i{:}], [args_o{:}], obj.class_name);

end

function str = get_signature(name,monom,sz)
    % Return signature.
    dim = get_dimensions(monom,sz);

    if ~isempty(dim)
        dim = {sprintf('[%s]',[dim{:}])};
    end

    str = sprintf('%s', name, [dim{:}]);
end

function dim = get_dimensions(monom,sz)
    % Return dimensions.
    n = sz(1); m = sz(2);
    d1 = monom.mindeg;
    d2 = monom.maxdeg;

    if n == 1 && m == 1
        % don't show dimensions
        dim = {};
    elseif m == 1
        % show number of rows
        dim = compose('%d',n);
    elseif n > 0 || m > 0
        % show size
        dim = compose('%dx%d',n,m);
    else
        % show empty dimensions
        dim = {''};
    end

    if all(d2 == 0)
        % nothing to do
        return
    elseif ~isempty(dim)
        dim{end+1} = ',';
    end

    if d1 == d2
        % show single degree
        dim{end+1} = sprintf('d=%d',d1);
    else
        % show degree range
        dim{end+1} = sprintf('d=%d:%d',d1,d2);
    end
end