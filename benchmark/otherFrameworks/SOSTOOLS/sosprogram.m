function sos = sosprogram(vartable,decvartable)
% SOSPROGRAM --- Initialize a new sum of squares program.
%
% SOSP = sosprogram(VARTABLE,DECVARTABLE)
%
% SOSP is a new sum of squares program.
% VARTABLE is a vector of independent variables.
% DECVARTABLE is a vector of decision variables (optional).
%
% Both VARTABLE and DECVARTABLE are either symbolic or polynomial objects.
%

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
% 12/24/01 - SP
% 01/07/02 - SP
% 02/21/02 - SP -- Symbolic polynomial
% 10/10/02 - SP -- Path checking
% 12/15/21 - DJ -- Correction for 'sos.var.idx{1}' specification in dpvar case
% 02/09/22 - DJ -- Add check for mosek solver
% 05/04/23 - DJ -- Convert 'decvartable' to row vector in 'sym' case.

if ~exist('sedumi') & ~exist('sqlp') & ~exist('csdp') & ~exist('sdpnal') & ~exist('sdpnalplus')&  ~exist('sdpam') &  ~exist('cdcs') & ~exist('mosekopt')
    error('No SDP solvers found.') ;
end;

if ~exist('getpolysym')
    dd = which('sosprogram');
    dd = strrep(dd,'/sosprogram.m','/');
    dd = strrep(dd,'\sosprogram.m','\');
    pp = genpath(dd);
    addpath(pp);
end;


if isa(vartable,'sym')
    
    sos.var.num = 0;
    sos.var.type = {};
    sos.var.Z = {};
    sos.var.ZZ = {};
    sos.var.T = {};
    
    sos.expr.num = 0;
    sos.expr.type = {};
    sos.expr.At = {};
    sos.expr.b = {};
    sos.expr.Z = {};
    
    sos.extravar.num = 0;
    sos.extravar.Z = {};
    sos.extravar.ZZ = {};
    sos.extravar.T = {};
    sos.extravar.idx = {};
    
    sos.objective = [];        % 01/07/02
    
    sos.solinfo.x = [];
    sos.solinfo.y = [];
    sos.solinfo.RRx = [];
    sos.solinfo.RRy = [];
    sos.solinfo.info = [];
    
    assume(vartable,'real');    % DJ - 05/04/23
    if ~isrow(vartable) %AP 30092020
        vartable = vartable.';
    end
	sos.vartable = converttochar(vartable);    % 30/09/2020 AP
    %     sos.vartable = sym2chartable(vartable);    % 02/21/02
%     if size(vartable,2) > 1
%         vartable = vartable.';
%     end;
    sos.symvartable = vartable;
    
    sos.varmat.vartable = '[]';                 %JA&GV 06/06/13
    sos.varmat.symvartable = [];
    sos.varmat.count = 0;
    
    
    if (nargin == 2 && ~isempty(decvartable))
        assume(decvartable,'real');     % DJ - 05/04/23
        if ~isrow(decvartable)
            decvartable = decvartable.';
        end

        sos.objective = sparse(length(decvartable),1);
        sos.decvartable = sym2chartable(decvartable);     % 03/01/02
%        setofCommaBrackets = [1 strfind(sos.decvartable,',') strfind(sos.decvartable,']')];%26/04/13
        
        if size(decvartable,2) > 1
            decvartable = decvartable.';%the transpose of a matrix of decision varibles
            %it is put under this form in
            %sos.symdecvartable (next line)
        end;
        sos.symdecvartable = decvartable;
        sos.var.idx{1} = length(find(sos.decvartable==','))+2;
        
    else
        sos.decvartable = '[]';
        sos.symdecvartable = [];
        sos.var.idx{1} = 1;
        
    end;
    
else
    
    sos.var.num = 0;
    sos.var.type = {};
    sos.var.Z = {};
    sos.var.ZZ = {};
    sos.var.T = {};
    
    sos.expr.num = 0;
    sos.expr.type = {};
    sos.expr.At = {};
    sos.expr.b = {};
    sos.expr.Z = {};
    
    sos.extravar.num = 0;
    sos.extravar.Z = {};
    sos.extravar.ZZ = {};
    sos.extravar.T = {};
    sos.extravar.idx = {};
    
    sos.objective = [];        % 01/07/02
    
    sos.solinfo.x = [];
    sos.solinfo.y = [];
    sos.solinfo.RRx = [];
    sos.solinfo.RRy = [];
    sos.solinfo.info = [];
    
    sos.vartable = sort(vartable.varname);
    
    sos.varmat.vartable = [];                 % PJS 9/9/2013
    sos.varmat.symvartable = [];
    sos.varmat.count = 0;
    
    
    if (nargin == 2 && ~isempty(decvartable))
        sos.objective = sparse(length(decvartable),1);
        if isa(decvartable,'dpvar')
            sos.decvartable = sort(decvartable.dvarname);
        else
            sos.decvartable = sort(decvartable.varname);
        end
	sos.var.idx{1} = length(sos.decvartable)+1;	% DJ, 12/15/21
    else
        sos.decvartable = {};
        sos.var.idx{1} = 1;
    end;
    
end;

