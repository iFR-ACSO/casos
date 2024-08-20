function [Q,Z,decomp,Den] = findsos(P,flag,options)
% FINDMATRIXSOS --- Find a sum of squares decomposition of a given matrix polynomial.
%
% [Q,Z,decomp,Den] = findsos(P,flag,options)
%
% P is a symmetric polynomial matrix.
%
% FLAG is an optional argument which, if set to 'rational', returns a
% sum of squares decomposition with rational coefficients. 
%
% OPTIONS is an optional argument which specifies the solver and the
% solver-specific parameters on its fields as
%   options.solver
%   options.params
% the default value for options.solver is 'SeDuMi'. The solver parameters
% are fields of options.params as, for instance, options.params.tol = 1e-9.
%
% A positive semidefinite Q and a symbolic monomial vector Z will be
% computed such that
%
%    (Ir kron Z)' * Q * (Ir kron Z) = P(x)
%
% If P is not a sum of squares, the function will return empty Q and Z.
%
% If P is a polynomial with integer coefficients and is represented as a
% symbolic object, then [Q,Z] = findsos(P,'rational') will compute a rational
% matrix Q such that Z'*Q*Z = P.

% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 4.00.
%
% Copyright (C)2002, 2004, 2013, 2016, 2018, 2021  
%                                      A. Papachristodoulou (1), J. Anderson (1),
%                                      G. Valmorbida (2), S. Prajna (3), 
%                                      P. Seiler (4), P. A. Parrilo (5),
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change log and developer notes


% 12/27/01 - SP
% 03/27/02 - SP
% 02/14/21 - DJ: Convert Q to full before taking sqrtm(Q)  in computing L
% 08/14/22 - DJ: Adjust polynomial implementation of 'rational' case. 
%                Also allow rational expansion in non-integer case.
% 09/03/22 - DJ: Initialize empty output fields if nargout>=3.
% 02/07/23 - DJ: Bugfix case where P has constant elements (on diagonal).
% 05/04/23 - DJ: Allow slightly negative eigenvalues in 'rational' case.

% Set tolerance for positivity of eigenvalues of Qr in rational case.
pos_tol_1 = -1e-14; % If min(eig(Qr)) <= pos_tol_1, return Qr with warning
pos_tol_2 = -1e-10; % If min(eig(Qr)) <= pos_tol_2, just return Q

if nargin == 2
    if ~strcmp(flag,'rational')
		options=flag;
		flag = 'abcdefgh';
	else
		options.solver='sedumi';
		%flag = 'abcdefgh';     % DJ, 08/14/22: should not overwrite flag, right?
    end
elseif nargin==1
	options.solver='sedumi';
end


% Check if matrix is square
dimp = size(P);
if dimp(1)~=dimp(2)
	disp('The polynomial matrix is not square, it cannot be a sum of squares');
	Q = [];
	Z = [];
    if nargout >= 3
		decomp = [];
		Den = [];
    end
	return
end

if isa(P,'double')
    % The object is just a constant matrix --> convert to polynomial
    P = polynomial(P);
end

if isa(P,'sym')
	
	P = expand(P);  % Do we need this? JA
	%vars = findsym(P);
	%vars = sym(['[',vars,']']);
	vars = symvar(P);  % JA Edit
	
	%     nvars = size(vars,2);
	%     for var = 1:nvars;
	%         for j = 1:dimp(1)
	%             if rem(double(feval(symengine,'degree',P(j,j))),2) ;
	%                 disp(['Degree in ' char(vars(var)) ' is not even for the diagonal elements. The polynomial matrix cannot be a sum of squares']);
	%                 Q = [];
	%                 Z = [];
	%                 return;
	%             end;
	%         end;
	%     end;
	%
	if isequal(P,P.')~=1;
		disp(['The polynomial matrix is not symmetric, it can not be a sum of squares']);
		Q = [];
		Z = [];
        if nargout >= 3
			decomp = [];
			Den = [];
        end
		return;
	end;
	
	prog = sosprogram(vars);
	prog = sosineq(prog,P);
    [prog,info] = sossolve(prog,options);
	
	if strcmp(options.solver,'SDPNAL')
		disp('findsos function currently not supported for SDPNAL solver');
		Q = [];
		Z = [];
		if nargout >= 3
			decomp = [];
			Den = [];
        end
		return;
	elseif (info.dinf==1)||(info.pinf==1)
		disp('No sum of squares decomposition is found.');
		Q = [];
		Z = [];
		if nargout >= 3
			decomp = [];
			Den = [];
        end
		return;
	else
		Q = reshape(prog.solinfo.RRx,sqrt(length(prog.solinfo.RRx)),sqrt(length(prog.solinfo.RRx)));
		if isa(P,'sym');
			vartable  = reshape(vars,1,length(vars));%makes vartable a row vector
			Z = prod(kron(vartable,ones(size(prog.extravar.Z{1},1),1)).^prog.extravar.Z{1},2);
		else
			pvar Z
			coefftemp = speye(size(prog.extravar.Z{1},1));
			Z = polynomial(coefftemp,prog.extravar.Z{1},prog.vartable,[size(prog.extravar.Z{1},1),1]);
		end;
	end;
	
	if nargin > 1 & strcmp(flag,'rational')
		
		A = (full(prog.expr.At{1})') ;
		B = (full(prog.expr.b{1})) ;
		% Round x0
		N = 1 ; % This is the tolerance in the rounding
		xx = prog.solinfo.x;
		
		kmax = 10 ;
		pd = 0 ;
		while kmax ~= 0 ;
			kmax = kmax - 1 ;
			x0 = round(N*xx);
            try [Qflat,NN] = proj3(x0,A,B,N);
                n = sqrt(length(Qflat));
                Qr = reshape(Qflat,n,n);
                % Qr should be PSD (should really check symbolically)
                if min(eig(Qr/NN))>=pos_tol_2 ; kmax=0 ; pd = 1 ; end
                % Increase N, and try again
                N = 2*N;
            catch
                NN = 1;
                n = sqrt(length(xx));
                Qr = reshape(xx,n,n);
                break
            end
        end
        % If eigenvalues of Qr are negative but small, return rational
        % decomposition with a warning.
        if pd==1 && min(eig(Qr/NN))<pos_tol_1    % DJ, 05/04/23
            fprintf(2,['Warning: Returned matrix Q in decomposition Zd''*Q*Zd has negative eigenvalue ',num2str(min(eig(Qr/NN))),'.\n',...
                       '         Rational decomposition may not be SOS.\n']);
        end
		
		% Experimental, no good error checking yet, so we check that
		% expand(NN*P - Z.'*Qr*Z) is zero!
		if (expand(NN*P-Z.'*Qr*Z) ~= 0) || (pd == 0)
% 			Qr = []; Z = []; NN = [];
% 			disp('Could not compute a rational SOS!');
            fprintf(2,['Warning: Could not compute a rational SOS decomposition;',...
                       ' returning real-valued decomposition instead.\n']);
            Qr = Q; NN = 1;
		end
		
		if nargout == 4
			Q = Qr;
			Den = NN;
		else
			if isa(P,'sym')
				Q = sym(Qr/NN);
			else
				Q = Qr/NN;
				disp('To obtain an exact rational solution, run the function with three output arguments.');
			end;
		end;
        
    elseif nargout==4
        Den = 1;
		
	end;
	
	
	L = real(sqrtm(full(double(Q)))); %AP edit for first order solver accuracy
	decomp = L*(kron(eye(dimp(1)),Z));
	
else 
	% PJS -- Handle polynomial variables
	% Most of this code simply mimics the symbolic case and hence the
	% overlap can be reduced.
	
	nvars = P.nvar;
	vars = polynomial(zeros(nvars,1));
	for j = 1:nvars
		vars(j) = pvar(P.varname{j});
	end
	
	for var = 1:nvars;
		for j = 1:dimp(1)
			Pjj = P(j,j);
            if isempty(Pjj.degmat)
                maxdeg = 0;
            else
			    maxdeg = max(Pjj.degmat(:,var));
            end
			if rem(maxdeg,2)
				disp(['Degree in ',vars(var).varname{1},' is not even for the diagonal elements. The polynomial matrix cannot be a sum of squares']);
				Q = [];
				Z = [];
                if nargout >= 3
                    decomp = [];
                    Den = [];
                end
				return;
			end;
		end;
	end;
	
	tmp = isequal(P,P.');
	if ~all(tmp(:))
		disp(['The polynomial matrix is not symmetric, it can not be a sum of squares']);
		Q = [];
		Z = [];
        if nargout >= 3
			decomp = [];
			Den = [];
        end
		return;
	end;
	
	prog = sosprogram(vars);
    if isempty(P.degmat)
        % DJ, 02/07/23: Add a decision var to the program to avoid error.
        [prog,~] = sossosvar(prog,1);
    end
    prog = sosineq(prog,dpvar(P));
	[prog,info] = sossolve(prog,options); %AP edit to pass solver.
	
	if strcmp('solver','SDPNAL')
		disp('findsos function currently not supported for SDPNAL solver');
		Q = [];
		Z = [];
		if nargout >= 3
			decomp = [];
			Den = [];
        end
		return;
	elseif (info.dinf==1)||(info.pinf==1)
		disp('No sum of squares decomposition is found.');
		Q = [];
		Z = [];
		if nargout >= 3
			decomp = [];
			Den = [];
        end
		return;
    else
        if isempty(P.degmat)
            % DJ, 02/07/23: Remove the added decision variable.
            Q = reshape(prog.solinfo.RRx(2:end),sqrt(length(prog.solinfo.RRx)-1),sqrt(length(prog.solinfo.RRx)-1));
        else
		    Q = reshape(prog.solinfo.RRx,sqrt(length(prog.solinfo.RRx)),sqrt(length(prog.solinfo.RRx)));
        end
		if isa(P,'sym')
			Z = mysympower(vars,prog.extravar.Z{1});
		else
			pvar Z
			coefftemp = speye(size(prog.extravar.Z{1},1));
			Z = polynomial(coefftemp,prog.extravar.Z{1},prog.vartable,[size(prog.extravar.Z{1},1),1]);
		end;
	end;
	
	if nargin > 1 & strcmp(flag,'rational')
		A = (full(prog.expr.At{1})') ;
		B = (full(prog.expr.b{1})) ;
		% Round x0
		N = 1 ; % This is the tolerance in the rounding
		xx = prog.solinfo.x;
		
		kmax = 10 ;
		pd = 0 ;
		while kmax ~= 0
			kmax = kmax - 1 ;
			x0 = round(N*xx);
            try [Qflat,NN] = proj3(x0,A,B,N);
                n = sqrt(length(Qflat));
                Qr = reshape(Qflat,n,n);
                % Qr should be PSD (should really check symbolically)
                if min(eig(Qr/NN))>=pos_tol_2 ; kmax=0 ; pd = 1 ; end
                % Increase N, and try again
                N = 2*N;
            catch
                NN = 1;
                n = sqrt(length(xx));
                Qr = reshape(xx,n,n);
                break
            end
        end

        % If eigenvalues of Qr are negative but small, return rational
        % decomposition with a warning.
        if pd==1 && min(eig(Qr/NN))<pos_tol_1    % DJ, 05/04/23
            fprintf(2,['Warning: Returned matrix Q in decomposition Zd''*Q*Zd has negative eigenvalue ',num2str(min(eig(Qr/NN))),'.\n',...
                       '         Rational decomposition may not be SOS.\n']);
        end
		
		% Experimental, no good error checking yet, so we check that
		% expand(NN*P - Z.'*Qr*Z) is zero!
		err = NN*P-Z.'*Qr*Z;
		if (isa(err,'polynomial') && max(max(abs(err.coeff)))>=1e-14) ...
                || (isa(err,'double') && ~all(err==0)) || (pd == 0)  % DJ, 08/14/22: remove call to "expand"
% 			Qr=[];Z=[];NN=[];
% 			disp('Could not compute a rational SOS!');
            fprintf(2,['Warning: Could not compute a rational SOS decomposition;',...
                       ' returning real-valued decomposition instead.\n']);
            Qr = Q; NN = 1;
		end
		
		if nargout == 4
			Q = Qr;
			Den = NN;
		else
			if isa(P,'sym')
				Q = sym(Qr/NN);
			else
				Q = Qr/NN;
				disp('To obtain an exact rational solution, run the function with three output arguments.');
			end;
		end;
        
    elseif nargout==4
        Den = 1;
		
	end;
		
	L = real(sqrtm(full(double(Q)))); %AP edit for first order solver accuracy
	decomp   = L*(kron(eye(dimp(1)),Z));
	
end;