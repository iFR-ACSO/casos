function [xs,ys,info] = sedumi2dsdp(As,bs,cs,K,opts) 
% function [x,y,info] = sedumi2dsdp(A,b,c,K,opts);
%
% DESCRIPTION 
%   This function solves a cone optimization specified in Sedumi format
%   using the DSDP solver.  The function converts the inputs from
%   Sedumi format to DSDP format, solves the problem with DSDP, and
%   converts the DSDP outputs back into Sedumi output format.
%
% INPUTS 
%   A,b,c,K:  Cone optimization problem in Sedumi format.  See 
%       Sedumi documentation for more details.  SOCP constraints
%       are not handled by DSDP.
%   opts: Options structure in DSDP format.

%       
% OUTPUTS
%   x,y: Optimization results in Sedumi output format. See Sedumi
%       documentation.
%   info: DSDP solver info.  See DSDP documentation.
%
% SYNTAX 
%   [x,y,info] = sedumi2dsdp(A,b,c,K,opts);
%
% EXAMPLE

% 2/28/08  PJS     Initial Coding

if ~isfield(K,'f')
    K.f = [];
end
if length(K.s)==1 && K.s==0
    K.s = [];
end

% DSDP does not handle free variables in the primal form. Split free 
% variabless into positive/negative parts.    
[Ad,bd,cd,Kd] = splitfreevars(As,bs,cs,K); 

% Convert from Sedumi Format to DSDP Format
bd = full(bd);
nblock = 1;
ptr = 0;
n = Kd.l;
if (n > 0)
    AC{nblock,1} = 'LP';
    AC{nblock,2} = n;
    AC{nblock,3} = [Ad(:,1:n)' cd(1:n)];
    nblock = nblock+1;
    ptr =n;
end
for i1 = 1:length(Kd.s)
    n = Kd.s(i1);
    AC{nblock,1} = 'SDP';
    AC{nblock,2} = n;
    idx = ptr+find(triu(reshape(1:n^2,n,n)));
    AC{nblock,3} = [Ad(:,idx)' cd(idx)];    
    nblock = nblock+1;
    ptr = ptr+n^2;
end

% Call DSDP (PDSDP?)

%if ~isfield(opts,'gaptol')
%    opts.gaptol = 1e-9;
%end
if ~isfield(opts,'print');
    opts.print = 0;
end
[info,yd,xd] = dsdp(bd,AC,opts);

% Convert solution from DSDP format to Sedumi format
ys = yd;
Ndl = (Kd.l>0);
if Ndl
    xs = xd{1};
else
    xs = [];
end
for i1=1:length(Kd.s);
    tmp = dmat(xd{i1+Ndl});
    xs = [xs; tmp(:)];
end

% Recombine positive/negative parts of the free variables 
if ~isempty(K.f)
    xs = [xs(K.f+1 : 2*K.f) - xs(1:K.f); xs(2*K.f+1:end)];    
end