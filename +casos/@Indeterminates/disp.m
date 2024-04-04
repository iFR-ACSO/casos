function disp(obj)
% Print indeterminate variables to command line.

s = str(obj);

if length(obj) == 1
    % single variable
    disp([s{:}])

    return
end

% combine to tuple
out = join(s,', ');

% display (x1,...,xN)
disp(['(' out{:} ')'])

end
