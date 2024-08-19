classdef Logger
% A simple logger class.

properties (Constant)
    % Predefined logger
    Error = casos.package.Logger(casos.package.VerbosityLevel.error,2);
    Debug = casos.package.Logger(casos.package.VerbosityLevel.debug,1);
    Off   = casos.package.Logger;
end

properties (Access=private)
    stream;
    vlevel;
end

methods
    function log = Logger(level,stream)
        % Create a new logger object.
        if nargin < 2
            % default stream: console
            stream = 1;
        end
        if nargin < 1
            % default verbosity: off
            level = casos.package.VerbosityLevel.off;
        end

        log.stream = stream;
        log.vlevel = casos.package.VerbosityLevel(level);
    end

    function varargout = printf(obj,lvl,fmt,varargin)
        % Write message to output stream if level permits.
        lvl = casos.package.VerbosityLevel(lvl);

        % check if requested level is within verbosity
        if lvl <= obj.vlevel
            % write to stream
            nb = fprintf(obj.stream,fmt,varargin{:});
        else
            % nothing to do
            nb = 0;
        end

        if nargout > 0
            varargout = {nb};
        end
    end
end

end
