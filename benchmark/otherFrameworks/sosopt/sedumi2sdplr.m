function [xs,ys,info] = sedumi2sdplr(As,bs,cs,K,opts)
% function [x,y,info] = sedumi2sdplr(A,b,c,K,opts);
%
% DESCRIPTION 
%   This function solves a cone optimization specified in Sedumi format
%   using the SDPLR solver.  The function converts the inputs from
%   Sedumi format to SDPLR format, solves the problem with SDPLR, and
%   converts the SDPLRoutputs back into Sedumi output format.
%
% INPUTS 
%   A,b,c,K:  Cone optimization problem in Sedumi format.  See 
%       Sedumi documentation for more details.  SOCP constraints
%       are not handled by SDPLR.
%   opts: Options structure in SDPLR format.
%       
% OUTPUTS
%   x,y: Optimization results in Sedumi output format. See Sedumi
%       documentation.
%   info: SDPLR solver info.  See SDPLR documentation.
%
% SYNTAX 
%   [x,y,info] = sedumi2sdplr(A,b,c,K,opts);
%
% EXAMPLE

% 11/15/10  PJS     Initial Coding

if ~isfield(K,'f')
    K.f = [];
end
if length(K.s)==1 && K.s==0
    K.s = [];
end

% Symmetrize equality constraints
As = symconstraints(As,K);

% SDPLR does not handle free variables in the primal form. Split free 
% variabless into positive/negative parts.  
[Ac,bc,cc,Kc] = splitfreevars(As,bs,cs,K); 

% Call SDPLR

%if ~isfield(opts,'objtol')
%    opts.feastol= 1e-9;
%end
if ~isfield(opts,'printlevel');
    opts.printlevel = 0;
end

% XXX PJS: On some probs SDPLR gives seg fault if data is sparse. This is 
% related to some corner-cases for sparse matrices in mexsdplr.c (max # of 
% nonzero entries ~= # of nonzero entries).  I think sparse A is handled 
% correctly and this only affects sparse b and c. For now convert b and
% c to full matrices.

%Ac = full(Ac);
bc = full(bc);
cc = full(cc);

[xs,ys,info,rr] = sdplr(Ac,bc,cc,Kc,opts);

% Recombine positive/negative parts of the free variables 
if ~isempty(K.f)
    xs = [xs(K.f+1 : 2*K.f) - xs(1:K.f); xs(2*K.f+1:end)];    
end

