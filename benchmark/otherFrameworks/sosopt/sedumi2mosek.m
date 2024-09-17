function [x, y, info] = sedumi2mosek(A, b, c, K, opts)
% function [x,y,info] = sedumi2mosek(A,b,c,K,opts);
%
% DESCRIPTION 
%   This function solves a cone optimization specified in Sedumi format
%   using the Mosek solver.  The function converts the inputs from
%   Sedumi format to Mosek format, solves the problem with Mosek, and
%   converts the Mosek outputs back into Sedumi output format.
%
% INPUTS 
%   A,b,c,K:  Cone optimization problem in Sedumi format.  See 
%       Sedumi documentation for more details.  
%   opts: Mosek solver options. Can include the fields 'cmd' and 'param' 
%         as described in the Mosek documentation.  Default for cmd is
%         'minimize echo(0)'.
%       
% OUTPUTS
%   x,y: Optimization results in Sedumi output format. See Sedumi
%       documentation.
%   info: Mosek solver info.  See Mosek documentation.
%
% SYNTAX 
%   [x,y,info] = sedumi2mosek(A,b,c,K,opts);


% Check options
mskcmd = [];
mskparam = [];
if nargin == 5 
    if isfield(opts, 'cmd')
        if  ~isempty(opts.cmd) && (~ischar(opts.cmd) || ~isrow(opts.cmd))
            error('cmd for mosekopt should be a charactar array');
        else
            mskcmd = opts.cmd;
        end
    end
    if isfield(opts, 'param') 
        if ~isstruct(opts.param)
            error('param for mosekopt should be a struct array');
        else
            mskparam = opts.param;
        end
    end
else
    error('Incorrect number of input arguments')
end
   
% Separate decision Variables
if isfield(K, 'f') % Free variables
    L(1) = K.f;
end
if isfield(K, 'l') % Nonegative variables
    L(2) = K.l;
end
if isfield(K, 'q') % Lorentz Cone
    error('Lorentz Cones are not currently implemented')
end
if isfield(K, 'r') % Rotated Lorentz Cone
    error('Lorentz Cones are not currently implemented')
end
if isfield(K, 's') % Semidefinite variable
    S = K.s;
    nS = cumsum([sum(L) S.^2]);
end

% Symmetrize equality constraints
A = symconstraints(A,K);

% Convert Problem
% Objective in linear vars
prob.c = c(1:sum(L));

% Objective in SD variables
Cmatmsk = [];
for j = 1:numel(S)
    Cmatmsk = mat2mosek(mat(c(1+nS(j):nS(j+1))), Cmatmsk, j);
end
prob.barc = Cmatmsk;

% Constraint Bounds
prob.blc = b;
prob.buc = b;

% Variable Bounds
prob.blx = -inf(sum(L),1);
prob.blx(L(1)+1:sum(L),1) = 0;
prob.bux = [];

% Constraint in scalar variables
prob.a = A(:,1:sum(L));

%Constraint in SD variables
prob.bardim = S; 
Amatmsk = [];
for j = 1:numel(S)
    for k = 1:size(A,1)
        Amatmsk = mat2mosek(mat(A(k,1+nS(j):nS(j+1))), Amatmsk, j, k);
    end
end
prob.bara = Amatmsk;

% Call Mosek
solveTime = tic;
[~,res] = mosekopt(['minimize echo(0) info ' mskcmd], prob, mskparam);  
solverTime = toc(solveTime);
xMos = res.sol.itr.barx;
x = zeros(sum(L)+nS(end),1);
x(1:sum(L)) = res.sol.itr.xx; 

% Reconstruct SD matrices to put in Sedumi format
initIdx = 1;
for k = 1:numel(S)
    xM = zeros(S(k));
    for j=1:S(k)
        xM(j:end,j) = xMos(initIdx:initIdx+S(k)-j);
        initIdx = initIdx+S(k)-j+1;
    end
    xM = xM + tril(xM,-1)'; 
    x(nS(k)+1:nS(k+1)) = vec(xM);
end
y = res.sol.itr.y;

% Create info structure
info = struct('pinf', 1, 'dinf', 1, 'mosekinfo', res, 'solvertime',solverTime );
switch res.sol.itr.solsta
    case {'OPTIMAL', 'NEAR_OPTIMAL', 'PRIM_AND_DUAL_FEAS', 'DUAL_FEAS', 'PRIM_FEAS'}
        info.pinf = 0;
        info.dinf = 0;
    case {'PRIM_INFEAS_CER'}
        info.dinf = 0;
    case {'DUAL_INFEAS_CER'}
        info.pinf = 0;
end

function mskmat = mat2mosek(matrix, mskmat, sdIdx, cnstIdx)
%  function mskmat = mat2mosek(matrix, mskmat, iIdx, jIdx)
%
%   Inputs:
%       matrix - 2-d matrix to add to mskmat structure in Mosek format
%       mskmat - stucture in the mosek sparse matrix format, if it is empty
%       it will be initialized
%           Fields: subi - index corresponding to constraint
%                   subj - index corresponding to semidefinite matrix
%                   subk - index corresponding to row
%                   subl - index corresponding to column
%                   val - value of entry
%       sdIdx - index of semidefinite matrix
%       cnstIdx[optional] - index of constraint  

if isempty(mskmat) && nargin == 4
    mskmat = struct('subi',[], 'subj', [], 'subk', [], 'subl', [], 'val', []);
elseif isempty(mskmat)
    mskmat = struct('subj', [], 'subk', [], 'subl', [], 'val', []);
end

[rIdx, cIdx, values] = find(tril(matrix));
if nargin == 4
    mskmat.subi = [mskmat.subi cnstIdx*ones(1, numel(values))];
end
mskmat.subj = [mskmat.subj sdIdx*ones(1, numel(values))];
mskmat.subk = [mskmat.subk rIdx'];
mskmat.subl = [mskmat.subl cIdx'];
mskmat.val = [mskmat.val values'];
