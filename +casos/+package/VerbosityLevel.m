classdef VerbosityLevel < uint32
% Level of verbosity for logging.

enumeration
    off     (0)     % no display
    error   (1)     % errors only
    warning (2)     % errors and warnings
    info    (3)     % information
    debug   (4)     % debug messages
    tout    (inf)   % everything
end

end
