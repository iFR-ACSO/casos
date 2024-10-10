classdef UnifiedReturnStatus

enumeration
    SOLVER_RET_SUCCESS
    SOLVER_RET_UNKNOWN
    SOLVER_RET_LIMITED
    SOLVER_RET_NAN
    SOLVER_RET_INFEASIBLE
end

methods
    function s = addfield(status,s)
        % Add unified return status to structure.
        s.UNIFIED_RETURN_STATUS = char(status);
    end
end

end
