% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [test_value,reference_value] = generateTestPolynomials(sz,variables,degmax,varargin)
% Generate test polynomials from CaΣoS and Multipoly

sp1 = generateRandomSparsity(sz,variables,degmax,varargin{:});

% generate random coefficients
coeffs1 = rand(1,sp1.nnz);
    
% build polynomials
poly1 = casos.PD(sp1,coeffs1);

% return test values as structs
arg1.sz = size(poly1);
[arg1.i,arg1.j,arg1.degrees] = get_tuplet(sparsity(poly1));
arg1.indets = str(poly1.indeterminates);
arg1.coeffs = coeffs1;

test_value = {arg1};

% convert to multipoly for reference
poly1_mltp = casos.toolboxes.to_multipoly(poly1);

% return reference values
reference_value = {poly1_mltp};

end
