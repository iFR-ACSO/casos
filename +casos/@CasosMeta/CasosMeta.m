classdef (Sealed) CasosMeta
% Meta information for CaΣoS.
%
% - version: Returns version in the format 'x.y.z(+)'.
% - git_revision: Returns git revision ID.
%

methods
    function obj = CasosMeta()
        error('Forbidden.')
    end
end

methods (Static)
    v = version;
    v = git_revision;
end

end