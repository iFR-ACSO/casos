%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% SPOTless. Implementation is based on the examples provided in the
% SPOTless toolbox and the manual
%
%--------------------------------------------------------------------------
function [gval,solverTime,buildTime]= reachEstGTM_benchSPOTless()

% indeterminates
x = msspoly ('x' , 4 ) ;
t = msspoly ('t' , 1 ) ;
x1 = x(1); % V 
x2 = x(2); % alpha
x3 = x(3); % q
x4 = x(4); % theta
% t = casos.Indeterminates('t');

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

solverTime1 = [];
solverTime2 = [];
% solverTime3 = [];
% bval_old = [];


% bisection tolerances (see default value of SOSOPT/GSOSOPT
relbistol = 1e-4;
absbistol = 1e-4;

   T = 3;
    h = t*(T-t);

%% V-s-iteration
% profile on -historysize 50000000000
for iter = 1:10

    % to make sure we do not use the old solution again
	gval = [];
    % bval = [];

    % solve gamma-step

    % find largest stable level set
    lb = 0; ub = 1;
 
    
    % bisection
    while (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
        % trial gamma
        gtry = (lb+ub)/2;
        
        % start time measure
        startTimeBuild1 = tic;
        
        % re-initialize a spot program
        pr1       = spotsosprog;
        pr1       = pr1.withIndeterminate(x);
        pr1       = pr1.withIndeterminate(t);

        % decision variable
        [pr1,s2]  = pr1.newFreePoly(monomials([x;t],0:4));
        [pr1,s3]  = pr1.newFreePoly(monomials([x;t],0:4));
        [pr1,s4]  = pr1.newFreePoly(monomials(x,0:4));
        [pr1,s5]  = pr1.newFreePoly(monomials([x;t],0:4));
        [pr1,s6]  = pr1.newFreePoly(monomials([x;t],0:4));
        [pr1,s7]  = pr1.newFreePoly(monomials([x;t],0:4));
        [pr1,s8]  = pr1.newFreePoly(monomials([x;t],0:4));
        [pr1,K]   = pr1.newFreePoly(monomials([x;t],0:4));
        


        % add sos constraints
        pr1 = pr1.withSOS( s2 );
        pr1 = pr1.withSOS( s3 );
        pr1 = pr1.withSOS( s4 - 0.0001);
        pr1 = pr1.withSOS( s5 );
        pr1 = pr1.withSOS( s6 );
        pr1 = pr1.withSOS( s7 );
        pr1 = pr1.withSOS( s8 );

        DVxval = diff(Vval, x);
        DVtval = diff(Vval, t);
        pr1 = pr1.withSOS(-(DVtval + DVxval*f + K*DVxval*g) - s2*h + s3*(Vval - gtry ) );
        pr1 = pr1.withSOS( -s4*rT + subs(Vval,t,T) - gtry  );
        pr1 = pr1.withSOS( uM - K + s5*(Vval - gtry ) - s6*h );
        pr1 = pr1.withSOS( K - um + s7*(Vval - gtry ) - s8*h );

         % solve problem
        opt         = spot_sdp_default_options();
        opt.verbose = 0;
        % opt.useQR = 1
        sol1 = pr1.minimize(msspoly(1),@spot_mosek,opt);
      
        % sol.status
        if strcmp(string(sol1.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE")
            % adapt lower interval bound
            lb = gtry;

            % store latest solution
            gval = gtry;
        
            s3val =  sol1.eval(s3);
            s5val =  sol1.eval(s5);
            s7val =  sol1.eval(s7);
            Kval  = sol1.eval(K);
        else
            % adapt upper interval bound
            ub = gtry;
        end

        
        % buildTime is total time spend to setup constraints (i.e sos problem),
        % do the transcription (poly --> sdp --> poly) we subtract the
        % solver time afterwards to only consider the actual build process
        endTimeBuild1 = [endTimeBuild1 toc(startTimeBuild1)-sol1.info.wtime];
        solverTime1   = [solverTime1 sol1.info.wtime];
		
    end

    sum( endTimeBuild1)
    sum(solverTime1)
    if ~isempty(gval)
        % fprintf('gamma is %g.\n', gval)
    else
         disp(['Problem infeasible in gamma-step in iteration:' num2str(iter)])
        break
    end
   
    

	% solve V-step
    % start time measure
	startTimeBuild2 = tic;
	
    % re-initialize sos program
    pr2       = spotsosprog;
    pr2       = pr2.withIndeterminate(x);
    pr2       = pr2.withIndeterminate(t);

    % decision variables
    [pr2,V]    =  pr2.newFreePoly(monomials([x;t],0:4));
    [pr2,s1]   =  pr2.newFreePoly(monomials(x,0:4));
    [pr2,s2]   =  pr2.newFreePoly(monomials([x;t],0:4));
    [pr2,s4]   =  pr2.newFreePoly(monomials(x,0:4));
    [pr2,s6]   =  pr2.newFreePoly(monomials([x;t],0:4));
    [pr2,s8]   =  pr2.newFreePoly(monomials([x;t],0:4));

    % constraints
    pr2 = pr2.withSOS( s1);
    pr2 = pr2.withSOS( s2);
    pr2 = pr2.withSOS( s6);
    pr2 = pr2.withSOS( s8);
    pr2 = pr2.withSOS( s4-0.0001  );
    DVx = diff(V, x);
    DVt = diff(V, t);
    pr2 = pr2.withSOS(-(DVt + DVx*f + Kval*DVx*g) - s2*h + s3val*(V - gval));
    pr2 = pr2.withSOS( -s4*rT + subs(V,t,T) - gval );
    pr2 = pr2.withSOS( uM - Kval + s5val*(V - gval) - s6*h );
    pr2 = pr2.withSOS( Kval - um + s7val*(V - gval) - s8*h );
    pr2 = pr2.withSOS( -(subs(V,t,0)-gval) + s1*(subs(Vval,t,0)-gval) );
    


    % solve problem
    opt         = spot_sdp_default_options();
    opt.verbose = 0;
    
    sol2 = pr2.minimize(msspoly(1),@spot_mosek,opt);
	
    % buildTime is total time spend to setup constraints (i.e sos problem),
    % do the transcription (poly --> sdp --> poly) we subtract the
    % solver time afterwards to only consider the actual build process
    endTimeBuild2(iter) = toc(startTimeBuild2)-sol2.info.wtime;
	solverTime2         = [solverTime2 sol2.info.wtime];
	
	if strcmp(string(sol2.info.solverInfo.itr.prosta),"PRIMAL_AND_DUAL_FEASIBLE")
           
            % extract solution
            Vval = sol2.eval(V);
	
			fprintf('Iteration %d: g = %g.\n',iter,full(gval));
	
	
	else
		disp(['Problem infeasible in V-step in iteration:' num2str(iter)])
		break
    end


	% check convergence
	% if ~isempty(bval_old)
    %     if abs(full(bval-bval_old)) <= 1e-3
    %         break
    %     else
    %         bval_old = bval;
    %     end
    % else
    %     bval_old = bval;
    % end


end % end for-loop
% profile viewer
buildTime  = sum(endTimeBuild1) + sum(endTimeBuild2);
solverTime = sum(solverTime1) + sum(solverTime2);
 

end % end of function
