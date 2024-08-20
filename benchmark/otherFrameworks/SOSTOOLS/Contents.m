% SOSTOOLS --- Sum of Squares Toolbox
% Version 4.00, 14 September 2021.
%
% Monomial vectors construction:
%    MONOMIALS   --- Construct a vector of monomials with prespecified 
%                    degrees.
%    MPMONOMIALS --- Construct a vector of multipartite monomials with
%                    prespecified degrees.
%
% General purpose sum of squares program (SOSP) solver:
% SOSPROGRAM        --- Initialize a new SOSP.
% SOSDECVAR         --- Declare new decision variables in an SOSP.
% SOSPOLYVAR        --- Declare a new polynomial variable in an SOSP.
% SOSSOSVAR         --- Declare a new sum of squares variable in an SOSP.
% SOSPOLYMATRIXVAR  --- Declare a new matrix of polynomial variables in an SOSP.
% SOSSOSMATRIXVAR   --- Declare a new matrix of sum of squares polynomial
%                       variables in an SOSP.
% SOSPOSMATR        --- Declare a new positive semidefinite matrix variable in an SOSP
% SOSPOSMATRVAR     --- Declare a new symbolic positive semidefinite matrix 
%                       variable in an SOSP
% SOSQUADVAR        --- Declare a polynomial/SOS decision variable with
%                       customized structure.
% SOSEQ             --- Add a new equality constraint to an SOSP.
% SOSINEQ           --- Add a new inequality constraint to an SOSP.
% SOSMATRIXINEQ     --- Add a new matrix inequality constraint to an SOSP.
% SOSSETOBJ         --- Set the objective function of an SOSP.
% SOSPSIMPLIFY      --- Perform simplification of an SOSP.
% SOSSOLVE          --- Solve an SOSP.
% SOSGETSOL         --- Get the solution from a solved SOSP.
%
% Customized functions:
%    FINDSOS     --- Find a sum of squares decomposition of a given polynomial.
%    FINDLYAP    --- Find a Lyapunov function for a dynamical system.
%    FINDBOUND   --- Find a global/constrained lower bound for a polynomial.
% 
% Demos:
%    SOSDEMO1 and SOSDEMO1S   --- Sum of squares test.
%    SOSDEMO2 and SOSDEMO2S   --- Lyapunov function search.
%    SOSDEMO3 and SOSDEMO3S   --- Bound on global extremum.
%    SOSDEMO4 and SOSDEMO4S   --- Matrix copositivity.
%    SOSDEMO5 and SOSDEMO5S   --- Upper bound for the structured singular value mu.
%    SOSDEMO6 and SOSDEMO6S   --- MAX CUT.
%    SOSDEMO7S                --- Chebyshev polynomials.
%    SOSDEMO8S                --- Bound in probability.
%    SOSDEMO9 and SOSDEMO9S   --- Matrix SOS decomposition.
%    SOSDEMO10 and SOSDEMO10S --- Set containment.

% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 4.00.
%
% Copyright (C)2002, 2004, 2013, 2016, 2018, 2021  A. Papachristodoulou (1), J. Anderson (1),
%                                      G. Valmorbida (2), S. Prajna (3), 
%                                      P. Seiler (4), P. A. Parrilo (5)
%                                      M. Peet (6), D. Jagt (6)
% (1) Department of Engineering Science, University of Oxford, Oxford, U.K.
% (2) Laboratoire de Signaux et Systmes, CentraleSupelec, Gif sur Yvette,
%     91192, France
% (3) Control and Dynamical Systems - California Institute of Technology,
%     Pasadena, CA 91125, USA.
% (4) Aerospace and Engineering Mechanics Department, University of
%     Minnesota, Minneapolis, MN 55455-0153, USA.
% (5) Laboratory for Information and Decision Systems, M.I.T.,
%     Massachusetts, MA 02139-4307
% (6) Cybernetic Systems and Controls Laboratory, Arizona State University,
%     Tempe, AZ 85287-6106, USA.
%
% Send bug reports and feedback to: sostools@cds.caltech.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
