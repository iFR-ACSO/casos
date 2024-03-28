function [A,B,ulim,xlim] = car2d(Ts)
% Model courtesy Ilya Kolmanovsky
m = 1891;	% mass
vx = 20;    % forward velocity
Cf = -18e3; % front tire cornering stiffness
Cr = -28e3; % rear tire cornering stiffness
lf = 1.5;   % front tire arm length
lr = 1.55;  % rear tire arm length
Iz = 3200;  % vehicle moment of inertia

% continuous-time system
A = [
        2*(Cf+Cr)/(m*vx)            2*(Cf*lf - Cr*lr)/(m*vx)-vx
        2*(Cf*lf - Cr*lr)/(Iz*vx)   2*(Cf*lf^2+Cr*lr^2)/(Iz*vx)
];
B = [-2*Cf/m; -2*Cf*lf/Iz];
C = eye(2);
D = 0;

xlim = [-3 3; -1 1];
ulim = [-0.15 0.15];

sysc = ss(A,B,C,D);

if nargin == 0 || Ts == 0
    % nothing to do
    return
end

% discrete-time system
sysd = c2d(sysc,Ts);
[A,B] = ssdata(sysd);
