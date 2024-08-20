function feas=checkfeas(info,sossol,opts)
% function feas=checkfeas(info,sossol,opts)
%
% DESCRIPTION
%   This function checks the feasibility of the solution of an SOS
%   optimization.  One or both of two feasibility checks can be performed:
%      1) Check the solver information contained in info for feasibility.
%      2) Check the Gram matrix decomposition in sossol for feasibility.
%
% INPUTS
%   info: Solution information structure returned by sossopt.
%   sossol: Solution cell array returned by sosopt.
%   opts: Options for optimization.  See SOSOPTIONS for more details.
%      Two options are used in the feasibility check:
%
%      -checkfeas: Check feasibility of solution.  'fast' checks
%         feasibility information returned by the solver. 'full' checks
%         the validity of the Gram matrix decomposition in the output
%         sossol. 'both' does both feasibility checks.
%         ['off'; {'fast'}; 'full'; 'both']
%      -feastol: Feasibility tolerance used in full feasibility check.
%         [{1e-6}]
%
% OUTPUTS
%   feas: feas=1 if the solution is feasible and 0 otherwise.
%
% SYNTAX
%   feas=checkfeas(info,sossol,opts)
%
% See also sosopt, sosoptions

% 4/25/2009 PJS Initial Coding
% 4/30/2009 PJS Added 'full' feasibility check
% 11/7/2010 PJS Pulled options logic into this file

% Parse inputs
if strcmpi(opts.solver,'setup')
    feas = 0;
    return;
elseif strcmpi(opts.checkfeas,'off');
    feas = 1;
    return;
elseif strcmpi(opts.checkfeas,'full');
    fastchk = 0;
    fullchk = 1;
elseif strcmpi(opts.checkfeas,'both');
    fastchk = 1;
    fullchk = 1;
elseif strcmpi(opts.checkfeas,'fast');
    fastchk = 1;
    fullchk = 0;
end
feastol = opts.feastol;

%-------- Feasibility check on solver information
fastfeas = 1;
if fastchk == 1
    % Get solver type and info returned by solver
    solverinfo = info.sdpsol.solverinfo;
    solver = opts.solver;
    
    % Check solution for feasibility
    if strcmpi(solver,'sedumi');
%         if  solverinfo.numerr ~= 2 && solverinfo.pinf == 0 && ...
%                 solverinfo.dinf==0 && solverinfo.feasratio > 0.70 ...
%                 && (solverinfo.numerr ~= 1 || solverinfo.feasratio>0.9)
        feasflg = (solverinfo.numerr == 0 && solverinfo.feasratio > 0.70) ...
            || (solverinfo.numerr==1 && solverinfo.feasratio> 0.9);

        if  solverinfo.pinf == 0 && solverinfo.dinf==0 && feasflg
            fastfeas = 1;
        else
            fastfeas = 0;
        end
    elseif strcmpi(solver,'dsdp')
        if strcmp(solverinfo.stype,'PDFeasible') && solverinfo.stopcode==1
            fastfeas = 1;
        else
            fastfeas = 0;
        end
    elseif strcmpi(solver,'sdplr');
        warning('NO FEAS CHECK FOR SDPLR----')
        fastfeas = 1;
    elseif strcmpi(solver,'sdpam');
        if strcmp(solverinfo.phasevalue,'pdOPT')
            fastfeas = 1;
        else
            fastfeas = 0;
        end
    elseif strcmpi(solver,'csdp');
        if solverinfo ==0 % || solverinfo==1 || solverinfo==2  || solverinfo ==3
            fastfeas = 1;
        else
            fastfeas = 0;
        end
    elseif strcmpi(solver,'sdpt3');
        % XXX What tolerances should be used here?
        if solverinfo.pinfeas<1 && solverinfo.dinfeas<1
            fastfeas = 1;
        else
            fastfeas = 0;
        end
    elseif strcmpi(solver,'mosek');
        if solverinfo.pinf<1 && solverinfo.dinf<1
            fastfeas = 1;
        else
            fastfeas = 0;
        end
    else
        fastfeas = 0;
    end
end

%-------- Feasibility check on Gram matrix decomposition
% XXX The full feas check needs further consideration.  For example,
% sosoptdemo5 currently fails the 'full' check because one of the error
% coefs is ce(206)= 3.3198e-008 and cs(206)=2.2653e-006.
% What is a good metric for the closeness of two polys?
% Should we have an absolute and relative tolerance?
fullfeas = 1;
if fullchk == 1
    % Check feasibility of solution
    Np = size(sossol,1);
    d = cell(Np,1);
    for i1=1:Np
        p = sossol(i1).p;
        z = sossol(i1).z;
        Q = sossol(i1).Q;
        
        if strcmp(p.RelOp,'==')
            % Check equality constraint
            s = p.OneSide;
            di = max(abs( s.coefficient ));
        else
            % Check (one-sided) SOS constraint
            s = p.OneSide;
            e = s - LOCALzQz(z,Q);
            
            % Project the coefs of s and e onto a common set of monomials
            M = monomials([s;e]);
            cs = abs(poly2basis(s,M));
            ce = abs(poly2basis(e,M));
            
            % Check equality: s=z'*Q*z
            idx1 = find(cs<=feastol);
            tmp1 = ce(idx1);
            
            idx2 = find(cs>feastol);
            tmp2 = ce(idx2)./cs(idx2);
            di(1) = max( [tmp1(:); tmp2(:)] );
            
            % Check Q >= 0
            reigQ=real(eig(Q));
            mine=min( reigQ );
            maxe=max( reigQ );
            di(2) = -mine/max(maxe,eps);
            
            
        end
        d{i1} = di(:);
    end
    
    d = cell2mat(d);
    if all(d<=feastol)
        fullfeas = 1;
    else
        fullfeas = 0;
    end
end

%-------- Total Feasibility
feas = fastfeas & fullfeas;

% if fastchk && fullchk && (fastfeas~=fullfeas)
%     warning('Disagreement in Feasibility -- Investigate the cause')
%     keyboard
% end



%-----------------------------------------------------------
% Local Function:
%   Fast implementation to perform  p = z'*Q*z
%   (This is now only slightly faster mtimes)
%-----------------------------------------------------------
function p = LOCALzQz(z,Q)

% Get z degmat with same ordering as the monomials appear in z
z = polynomial(z(:));
zdeg = z.coefficient'*z.degmat;
lz = length(z);

% Create a lower triangular Qlt such that z'*Qlt*z = z'*Q*z
% [This will allow us to express the degmat with half the number of terms]
Qlt = tril(Q)+triu(Q,+1)';

% Create coefficient vector
pcoef = Qlt(:);
idx = find(Qlt~=0);
pcoef = pcoef(idx);

if isempty(pcoef)
    % PJS 1/19/2012: All entries of Q are zero-->p=0
    p = polynomial(0);
else
    % Compute degmat for p = z*Q*z
    nxvar = size(zdeg,2);
    pdeg = repmat(zdeg,[lz 1]) +reshape(repmat(zdeg(:)',[lz 1]),[lz^2 nxvar]);
    pdeg = pdeg(idx,:);

    % Create polynomial
    pvarname = z.var;
    p = polynomial(pcoef,pdeg,pvarname,[1 1]);
    p = combine(p);
end