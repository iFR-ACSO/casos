% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function res = eval_vertcat(arg1,arg2)
% Evaluate vertical concatenation.

% detect empty polynomials
if (isempty(arg1) && isempty(arg2))
    if (size(arg1,2) == 0 && size(arg2,2) == 0)
        sz = size(arg1) + size(arg2);
    else
        sz = [
            min(size(arg1,1),size(arg2,1))
            max(size(arg1,2),size(arg2,2))
        ]';
    end
    
    res = polynomial(zeros(sz));
    return

elseif isempty(arg1)
    res = arg2;
    return

elseif isempty(arg2)
    res = arg1;
    return
end

% else
res = vertcat(arg1,arg2);

end
