% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function ref = multipoly2struct(poly)
% Convert a multipoly reference polynomial to a structure.

sz = size(poly);

if isempty(poly)
    % convert empty polynomial
    ref = struct('sz',sz,'i',zeros(1,0),'j',zeros(1,0));
    ref.degrees = [];
    ref.indets = {};
    ref.coeffs = zeros(0,1);
    return

else
    % ensure polynomial without zero terms
    poly = polynomial(poly + 0);
end

coeffs = poly.coefficient;
degmat = poly.degmat;

% sort monomials to graded REVERSE lexicographic order
% see casos.Sparsity/uniqueDeg
nt = size(coeffs,1);
ne = size(coeffs,2);

% sort by ascending degree
degsum = sum(degmat,2);

% make degree matrix unique
[degmat2,id,ic] = unique([degsum fliplr(degmat)],'rows','sorted');

% reverse order of degrees
degmat = fliplr(degmat2(:,2:end));

% sum repeated coefficients
summat = sparse(ic,1:nt,1,length(id),nt);
coeffs = (summat*coeffs);

% get tuplet
[row,col] = find(coeffs);
[i1,j1] = ind2sub(sz,reshape(col,1,[])); % ensure row vectors

% create struct
ref = struct('sz',sz,'i',i1-1,'j',j1-1); % 0-index

% store variables and degrees
ref.degrees = full(degmat(row,:));
ref.indets = poly.varname;
% store nonzero coefficients
ref.coeffs = nonzeros(coeffs);

end
