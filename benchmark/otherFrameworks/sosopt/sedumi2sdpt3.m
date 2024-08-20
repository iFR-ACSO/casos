function [xs,ys,info] = sedumi2sdpt3(As,bs,cs,K,opts)
% function [x,y,info] = sedumi2sdpt3(A,b,c,K,opts);
%
% DESCRIPTION 
%   This function solves a cone optimization specified in Sedumi format
%   using the SDPT3 solver.  The function converts the inputs from
%   Sedumi format to SDPT3 format, solves the problem with SDPT3, and
%   converts the SDPT3 outputs back into Sedumi output format.
%
% INPUTS 
%   A,b,c,K:  Cone optimization problem in Sedumi format.  See 
%       Sedumi documentation for more details.  
%   opts: Options structure in SDPT3 format.
%       
% OUTPUTS
%   x,y: Optimization results in Sedumi output format. See Sedumi
%       documentation.
%   info: SDPT3 solver info.  See SDPT3 documentation.
%
% SYNTAX 
%   [x,y,info] = sedumi2sdpt3(A,b,c,K,opts);
%
% EXAMPLE

% 9/18/09  PJS     Initial Coding

if ~isfield(K,'f')
    K.f = [];
end
if length(K.s)==1 && K.s==0
    K.s = [];
end

% Options
if ~isfield(opts,'printlevel');
% OPTIONS.printlevel  : 3, if want to display result in each iteration, 
%                       2, if want to display only summary,
%                       1, if want to display warning message,
%                       0, no display at all.  
    opts.printlevel = 0;
end

% Symmetrize equality constraints: 
% SQLP/read_sedumi will symmetrize but it also prints a warning
As = symconstraints(As,K);

% Convert to SDPT3 format
[blk,At,C,b] = read_sedumi(As,bs,cs,K);

% Solve with SDPT3
[obj,X,y,Z,info] = sqlp(blk,At,C,b,opts); 

% Convert solution back to Sedumi format
xs = [];
for i1=1:length(X)
    if strcmp( blk{i1,1}, 's') 
        Ns = length( blk{i1,2} );
        ptr = [0 cumsum(blk{i1,2})];
        for i2=1:Ns
            idx = ptr(i2)+1:ptr(i2+1);
            Xblk = X{i1}(idx,idx);             
            xs = [xs;Xblk(:)];
        end
    else
        xs = [xs;X{i1}(:)];
    end
end
xs = full(xs);
ys = full(y);

