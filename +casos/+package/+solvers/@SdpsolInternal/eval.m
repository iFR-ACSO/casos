function argout = eval(obj,argin)
% Call SDP interface.

data = call(obj.fhan,argin);

% to double
% data = cellfun(@sparse, data, 'UniformOutput',false);

% call conic solver
sol = call(obj.solver,data);

% parse solution
argout = call(obj.ghan,sol);

end
