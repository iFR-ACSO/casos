clc
clear

%% define variables
x = casos.PS('x',6,1);
u = casos.PS('u',3,1);


x1 = x(1); 
x2 = x(2); 
x3 = x(3); 
x4 = x(4);
x5 = x(5);
x6 = x(6);

J = diag([8970;9230;3830]);% Casini Parameter

% skew-symmetric matrix
skew = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];


B = @(sigma) (1-sigma'*sigma)*eye(3)+skew(sigma)+ 2*sigma*sigma';

% dynamics
x_dot =  [-inv(J)*skew(x(1:3))*J*x(1:3) + inv(J)*u; % omega_dot
           1/4*B(x(4:6))*x(1:3)];                     % MRP kinematics

% just for testing set to identity matrices
Dx = eye(6);
Du = eye(3);

f_scaled = Dx*subs(x_dot,[x;u],[Dx^-1*x;Du^-1*u])
