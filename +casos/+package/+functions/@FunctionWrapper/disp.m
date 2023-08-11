function disp(obj)
% Display function signature.

args_i = cellfun(@get_signature, struct2cell(obj.wrap.arg_i), 'UniformOutput', false);
args_o = cellfun(@get_signature, struct2cell(obj.wrap.arg_o), 'UniformOutput', false);

% separate arguments by comma
args_i(2,1:end-1) = {','};
args_o(2,1:end-1) = {','};

fprintf('%s:(%s)->(%s) %s\n', obj.name, [args_i{:}], [args_o{:}], obj.class_name);

end