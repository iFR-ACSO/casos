function disp(obj)
% Display function signature.

args_i = cell(1,obj.n_in);
args_o = cell(1,obj.n_out);

for idx=0:obj.n_in-1
    args_i{idx+1} = get_signature(name_in(obj,idx),sparsity_in(obj,idx));
end
for idx=0:obj.n_out-1
    args_o{idx+1} = get_signature(name_out(obj,idx),sparsity_out(obj,idx));
end

% separate arguments by comma
args_i = join(args_i,',');
args_o = join(args_o,',');

fprintf('%s:(%s)->(%s) %s\n', obj.name, args_i{:}, args_o{:}, obj.class_name);

end

function str = get_signature(name,sp)
    % Return signature.
    dim = signature(sp);

    if ~isempty(dim)
        dim = compose('[%s]',[dim{:}]);
    end

    str = sprintf('%s', name, [dim{:}]);
end
