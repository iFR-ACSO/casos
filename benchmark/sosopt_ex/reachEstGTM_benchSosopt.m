%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% SOSOPT/GSOSOPT. GSOSOPT performs the bisection.
% Implementation is based on the example gsosoptdemo1.m provide in the
% SOSOPT toolbox.
%
%--------------------------------------------------------------------------
function [gval,solverTime,buildTime,Vval]= reachEstGTM_benchSosopt()

% System dynamics
x = mpvar('x',[4,1]);
x1 = x(1); % V 
x2 = x(2); % alpha
x3 = x(3); % q
x4 = x(4); % theta

pvar t b

% Polynomial Dynamics
d2r = pi/180;
% scale matrix
Dmax = diag([20 20*d2r 50*d2r 20*d2r]);

% scaled 4-state system
f = [ % f1
      - 0.01955829859645207*x1^2*x2 + 0.0006116050110168922*x1^2*x3 - 0.4597074905072323*x1*x2^2 ...
      - 0.02143363124979007*x1*x2*x3 + 0.0913633506074555*x2^3 + 0.0104276150041391*x2^2*x4 ...
      - 0.0104276150041391*x2*x4^2 + 0.003475871668046367*x4^3 - 0.00975287614756193*x1^2 ...
      - 0.08801234368403432*x1*x2 - 0.001251382140034897*x1*x3 - 0.5201826057909725*x2^2 ...
      - 0.04982793763401992*x2*x3 + 0.00180695498489715*x3^2 - 0.04388794266402869*x1 ...
      + 0.07194181873717953*x2 - 0.00290915666614335*x3 - 0.1711592037553279*x4;...
      % f2
      0.05879789171572614*x1^3 + 0.6755306062430397*x1^2*x2 + 0.07878176650294093*x1^2*x3 ...
      + 0.6715746030304725*x1*x2^2 - 0.03631259653379534*x1*x2*x4 - 0.006619604539373272*x1*x3^2 ...
      + 0.01815629826689767*x1*x4^2 - 0.1817374041405456*x2^3  + 0.1364168103441202*x1^2 ...
      - 1.371702883566865*x1*x2 + 0.005525130555690897*x1*x3 + 1.479555253341581*x2^2 ...
      + 0.002672981895904953*x2*x3 + 0.07915793925313672*x2*x4 + 0.0135837161533212*x3^2 ...
      - 0.03957896962656836*x4^2  - 0.5767816606949114*x1 - 3.236059303865032*x2 + 2.30669948822443*x3;...
      % f3
      - 3.582346855286648*x1^2*x2 + 0.9194217971573437*x1^2*x3 + 2.279359758657476*x1*x2^2 ...
      - 0.3522980647489417*x2^3 - 0.07456233935706458*x1^2 - 16.12056084878992*x1*x2 ...
      - 1.881194554322757*x1*x3 + 2.56427972848966*x2^2  - 0.33553052710679*x1 ...
      - 18.13563095488866*x2 - 4.373316114187153*x3;...
      % f4
      2.5*x3];
g = [% g1
     - 0.01200284890938563*x1^2- 0.0672488718061946*x1*x2- 0.05401282009223534*x1 ...
     - 0.07565498078196893*x2- 0.06076442260376476; ...
     % g2
     0.190652302942119*x1^2 + 0.1033082169845372*x1*x2 - 0.3900865469718997*x1 ...
     + 0.240166275746567*x2 - 0.9068555816726751; ...
     % g3
     -13.57875281932966*x1^2 + 14.75731515054161*x1*x2 - 61.10438768698349*x1 ...
     + 16.60197954435931*x2 - 68.74243614785642; ...
     % g4
     0
    ];

% target function
rT = x'*Dmax'*blkdiag(1/(4)^2, 1/(pi/30)^2, 1/(pi/15)^2, 1/(pi/30)^2)*Dmax*x - 1;


P = [4.54064415441163	0.139789072082587	-0.0107315729582360	-0.406387155887787
    0.139789072082587	0.165595957230149	-0.00115647767230559	-0.0841739454608939
-0.0107315729582360	-0.00115647767230559	0.00732884319609208	0.00922486739143782
-0.406387155887787	-0.0841739454608939	0.00922486739143782	0.514580693591350];

Vval = x'*P*x;

% Trim point for elevator channel is 0.0489 rad.
% Saturation limit for elevator channel is -10 deg to 10 deg
uM = 10*d2r - 0.0489;
um = -(10*d2r + 0.0489);



% initialize arrays
endTimeBuild1 = [];
endTimeBuild2 = [];
gval_old = [];

solverTime1 = zeros(100,1);
solverTime2 = zeros(100,1);


% decision variables
K  = polydecvar('k',monomials([x;t], 0:4 ));
V  = polydecvar('v',monomials([x;t], 0:4 ));
s1 = sosdecvar('s1', monomials(x,0:2));
s2 = sosdecvar('s2',monomials([x;t],0:2));
s3 = sosdecvar('s3', monomials([x;t],0:2));
s4 = sosdecvar('s4', monomials(x,0:2));
s5 = sosdecvar('s5', monomials([x;t],0:2));
s6 = sosdecvar('s6', monomials([x;t],0:2));
s7 = sosdecvar('s7', monomials([x;t],0:2));
s8 = sosdecvar('s8', monomials([x;t],0:2));

T = 3;
h = t*(T-t);

DVx = jacobian(V, x);
DVt = jacobian(V, t);


%% V-s-iteration
for iter = 1:10
  
    opts        = gsosoptions;
    opts.minobj = -1; 
    opts.maxobj = 0;
    opts.absbistol = 1e-4;
    opts.relbistol = 1e-4;
    opts.solver = 'mosek';
    % opts.simplify = 'off';
  

    % gamma-step    
    startTimeBuild1 = tic;

    DVxval = jacobian(Vval, x);
    DVtval = jacobian(Vval, t);

    sosc    = [s2;s3;s4-0.0001;s5;s6;s7;s8;...
            -(DVtval + DVxval*f + K*DVxval*g) - s2*h + s3*(Vval + b); ...
            -s4*rT + subs(Vval,t,T) + b;...
            uM - K + s5*(Vval + b) - s6*h; ...
            K - um + s7*(Vval + b) - s8*h];
    
    % Solve with gsosopt
    [info,dopt] = gsosopt(sosc,[x;t],b,opts);

     % extract solution
      Kval = subs(K,dopt);
      s3val = subs(s3,dopt);
      s5val = subs(s5,dopt);
      s7val = subs(s7,dopt);
      gval  = -info.tbnds(2);
    
    % time measure of first subproblem
    solverTime1(iter) = sum(arrayfun(@(x) x.solverinfo.solvertime, info.sdpsol));

    % buildTime is total time spend to setup constraints (i.e sos problem),
    % do the transcription (poly --> sdp --> poly); in gsosopt the problem
    % is automatically re-build in each subiteration; we subtract the
    % solver time afterwards to only consider the actual build process
    endTimeBuild1     = [endTimeBuild1 toc(startTimeBuild1)-solverTime1(iter)]; 


    % V-step    
    startTimeBuild2 = tic;

    sosc    = [s1;s2;s4-0.0001;s6;s8;...
            -(subs(V,t,0)-gval) + s1*(subs(Vval,t,0)-gval);
            -(DVt + DVx*f + Kval*DVx*g) - s2*h + s3val*(V - gval);...
            -s4*rT + subs(V,t,T) - gval;...
            uM - Kval + s5val*(V - gval) - s6*h; ...
            Kval - um + s7val*(V - gval) - s8*h];
    
    % Solve with sosopt
    opts = sosoptions;
    opts.solver = 'mosek';
    % opts.simplify = 'off';
      
    [info,dopt,~] = sosopt(sosc,[x;t],opts);

    % time measure of third subproblem
    solverTime2(iter)   = info.sdpsol.solverinfo.solvertime;

    % buildTime is total time spend to setup constraints (i.e sos problem),
    % do the transcription (poly --> sdp --> poly); we subtract the
    % solver time afterwards to only consider the actual build process
    endTimeBuild2(iter) = toc(startTimeBuild2)-solverTime2(iter);
    
    % extract solution
    Vval = subs(V,dopt);

    fprintf('Iteration %d: g = %g.\n',iter ,full(gval));
	
				
	% check convergence
	% if ~isempty(gval_old)
		% 	if abs(full(gval-gval_old)) <= 1e-3
			% 	break
		% 	else
			% 	gval_old = gval;
		% 	end
	% 	else
		% 	gval_old = gval;
	% end

end % end-for-loop


buildTime  = sum(endTimeBuild1) + sum(endTimeBuild2);
solverTime = sum(solverTime1)   + sum(solverTime2) ;


% save the complete workspace, so people do not have to re-run execpt they
% want to
% save('SOSOPT_GTM_reach_bench.mat')

end % end of function