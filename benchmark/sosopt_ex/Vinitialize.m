function P = Vinitialize
syms x1 x2 x3 x4 u
x = [x1; x2; x3; x4];
nX = size(x,1);
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
xdot = f + g*u;
Apre = jacobian(xdot,x);
Bpre = jacobian(xdot,u);

% equilibriums
xbar = zeros(nX,1);
ubar = 0;

% substitute in the value of equilibrium
A = double(subs(subs(Apre,x,xbar),u,ubar));
B = double(subs(subs(Bpre,x,xbar),u,ubar));

% Q = diag([20 10 20 10 20 10]);
% R = diag([100 100]);
% K = lqr(A,B,Q,R);
% 
% P = lyap((A-B*K)',10*eye(6));

Q = diag([1 1 1 1]);
R = diag(1);
K = lqr(A,B,Q,R);

P = lyap((A-B*K)',eye(nX));