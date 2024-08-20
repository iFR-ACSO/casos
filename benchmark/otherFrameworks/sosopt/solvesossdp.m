function [x,y,solverinfo] = solvesossdp(sdpdata,opts);
% function [x,y,solverinfo] = solvesossdp(sdpdata,opts);
%
% DESCRIPTION
%   This function solves a cone optimization specified in Sedumi format
%   using the specified by solver.  Solver can be 'dsdp', 'csdp'
%   'sdpam', 'sedumi', 'setup' (Default is 'sedumi')
%
% INPUTS
%   sdp: Structure describing the cone optimization problem in Sedumi
%       format.  The fields are A,b,c,K.  See Sedumi
%       documentation for more details on the format.
%   opts: Structure of options with the following fields:
%       -solver: Optimization solver to be used. Choices are 'sedumi',
%          'sdpam', 'csdp', 'dsdp', or 'setup'.  solver='setup' will
%          simply return empty solutions. (Default is 'sedumi'.)
%       -The remaining fields of opts are passed directly to the solver.
%         Default is to set the solver options to turn off display.
%
% OUTPUTS
%   x,y: Optimization results in Sedumi output format. See Sedumi
%       documentation.
%   solverinfo: Info returned by solver.  See documentation for the
%       specified solver.
%
% SYNTAX
%   [x,y,solverinfo] = solvesossdp(sdpdata);
%   [x,y,solverinfo] = solvesossdp(sdpdata,opts);
%
% EXAMPLE
%

% 5/5/09  PJS  Initial Coding

% Get sdp problem data
A=sdpdata.A;
b=sdpdata.b;
c=sdpdata.c;
K=sdpdata.K;

% Set default options
if nargin<2
    opts = sosoptions;
elseif ~isa(opts,'sosoptions')
    error('See SOSOPTIONS for creation of options object.');
end
solver = opts.solver;
solveropts = opts.solveropts;

% Solve SDP
if exist(['sedumi2' solver], 'file')
    % call requested solver if conversion from Sedumi is supported
    
    % XXX Sedumi allows for non-symmetric constraints, e.g. introduce
    % a 2-by-2 Q with a PSD constraint and then add an equality constraint
    % Q(1,2) = 4.  SDPT3 gives a warning and converts this to a
    % symmetric constraint.
    [x,y,solverinfo] = feval(['sedumi2' solver],A,b,c,K,solveropts);
    return
end

% else:
switch (solver)
case 'sedumi'
    %if ~isfield(solveropts,'eps')
    %    solveropts.eps = 1e-9;
    %end
    if ~isfield(solveropts,'fid')
        solveropts.fid = 0;
    end
        
    if length(c)<= length(b)
        % 11/11/09 -- Ufuk found that sosopt(x^2) causes an error. This
        % is because sedumi errors if length(c) <= length(b). Add dummy
        % variables for now to catch this case.
        Npad = length(b)-length(c)+1;
        cpad = [c; zeros(Npad,1)];
        Apad = [A  zeros(size(A,1),Npad)];
        Kpad = K;
        Kpad.s = [Kpad.s ones(1,Npad)];
        [xpad,y,solverinfo] = sedumi(Apad,b,cpad,Kpad,solveropts);
        x = xpad(1:end-Npad);
%     elseif isempty(b)
%         % 11/29/10 -- sosopt(x^2) in kernel form gives an SDP with
%         % both A and b as 0-by-1. This causes Sedumi to error out.
%         % sosopt(x1^2+x2^2) and sosopt(5) also cause errors.
%         b = 0;
%         Npad = length(b)-length(c)+1;
%         cpad = [c; zeros(Npad,1)];
%         Apad = [A  zeros(size(A,1),Npad)];
%         Kpad = K;
%         Kpad.s = [Kpad.s ones(1,Npad)];
%         [xpad,y,solverinfo] = sedumi(Apad,b,cpad,Kpad,solveropts);
%         x = xpad(1:end-Npad);
%         y = [];        
    else
        [x,y,solverinfo] = sedumi(A,b,c,K,solveropts);
    end
    if isempty(x)
        % 11/15/10 -- LC found an example with infeasible eq. constraints.
        % Sedumi returns x as empty and solverinfo only contains the
        % field solver.pinf= 1.
        x = zeros(size(A,2),1);
        solverinfo.iter = 0;
        solverinfo.feasratio = -1;
        solverinfo.pinf = 1;
        solverinfo.dinf = 1;
        solverinfo.numerr = 0;
        solverinfo.timing = 0;
        solverinfo.cpusec = 0;
    end
    
case 'setup'
    x=[]; y=[]; solverinfo=[];
    
otherwise
    error(['Solver ' solver ' is not available']);
end

end