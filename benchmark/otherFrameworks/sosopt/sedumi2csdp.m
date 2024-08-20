function [xs,ys,info] = sedumi2csdp(As,bs,cs,K,opts)
% function [x,y,info] = sedumi2csdp(A,b,c,K,opts);
%
% DESCRIPTION 
%   This function solves a cone optimization specified in Sedumi format
%   using the CSDP solver.  The function converts the inputs from
%   Sedumi format to CSDP format, solves the problem with CSDP, and
%   converts the CSDP outputs back into Sedumi output format.
%
% INPUTS 
%   A,b,c,K:  Cone optimization problem in Sedumi format.  See 
%       Sedumi documentation for more details.  SOCP constraints
%       are not handled by CSDP.
%   opts: Options structure in CSDP format.
%       
% OUTPUTS
%   x,y: Optimization results in Sedumi output format. See Sedumi
%       documentation.
%   info: CSDP solver info.  See CSDP documentation.
%
% SYNTAX 
%   [x,y,info] = sedumi2csdp(A,b,c,K,opts);
%
% EXAMPLE

% 2/28/08  PJS     Initial Coding

if ~isfield(K,'f')
    K.f = [];
end
if length(K.s)==1 && K.s==0
    K.s = [];
end

% CSDP does not handle free variables in the primal form. Split free 
% variabless into positive/negative parts.  
% XXX Use CSDP function CONVERTF to perform this variable splitting?
[Ac,bc,cc,Kc] = splitfreevars(As,bs,cs,K); 

% Call CSDP

%if ~isfield(opts,'objtol')
%    opts.objtol= 1e-9;
%end
if ~isfield(opts,'printlevel');
    % XXX CSDP still prints out four lines of basic timing info even
    % if printlevel = 0. Can this be turned off?
    opts.printlevel = 0;
end
Ac = full(Ac);
bc = full(bc);
cc = full(cc);
[xs,ys,zs,info] = csdp(Ac,bc,cc,Kc,opts);

% Recombine positive/negative parts of the free variables 
if ~isempty(K.f)
    xs = [xs(K.f+1 : 2*K.f) - xs(1:K.f); xs(2*K.f+1:end)];    
end

