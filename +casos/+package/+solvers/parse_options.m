function [options,opts] = parse_options(options,opts)
% Parse options.
    if nargin < 2
        opts = struct;
    end
    for fn=intersect(fieldnames(options),fieldnames(opts))'
        if isempty(fn), continue; end
        options.(fn{:}) = opts.(fn{:});
        opts = rmfield(opts,fn);
    end
end
