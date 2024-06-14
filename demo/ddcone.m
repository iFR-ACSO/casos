% Test and verify DD cone 

% dimension of SDP
ndim = 5;

% decision variables
x = casadi.SX.sym('x', 2, 1);
X = casadi.SX.sym('X', ndim, ndim);

% create random SDP data
A = randn(ndim, ndim);
B = randn(ndim, ndim);
A = A+A';
B = B+B';

% define SDP problem
sdp.f = x(1);
sdp.g = eye(ndim) + x(1)*A + x(2)*B - X;
sdp.g = sdp.g(:);
sdp.x = [x; X(:)];

% define cones
opts = struct;
opts.Kx = struct('lin', 2, 'psd', ndim);
opts.Kc = struct('lin', length(sdp.g));

% initialize solver
S = casos.sdpsol('S','sedumi',sdp,opts);

% solve with equality constraints
tic
sol = S('lbg', zeros(ndim^2,1),'ubg',zeros(ndim^2,1),'lbx',-inf,'ubx',+inf);
fprintf('time with psd: %d \n', toc);

disp(S.stats.UNIFIED_RETURN_STATUS);
fprintf('solution with psd: %d \n', sol.f.full)

% Solve with DD in the decision variables
clear opts
opts.Kx = struct('lin', 2, 'dd', ndim);
opts.Kc = struct('lin', length(sdp.g));
% initialize solver
S = casos.sdpsol('S','sedumi',sdp,opts);

% solve with equality constraints
tic
sol = S('lbg', zeros(ndim^2,1),'ubg',zeros(ndim^2,1),'lbx',-inf,'ubx',+inf);
fprintf('time with dd in the decision variables: %d \n', toc);
disp(S.stats.UNIFIED_RETURN_STATUS);
fprintf('solution with dd: %d \n', sol.f.full)


%% plot the set of feasible solutions
xp = -0.5:0.01:0.5;
yp = -0.5:0.01:0.5;
[xp, yp] = meshgrid(xp,yp);
pplot = [];
dplot = [];
for i=1:size(xp,1)
    for j=1:size(xp,2)
        % verify if it is psd
        deig = eig(eye(ndim) + xp(i,j)*A + yp(i,j)*B);
        if all(deig>=0)
            pplot = [pplot; xp(i,j) yp(i,j)];
        end
        % verify if it is diagonally dominant
        if all(2*abs(diag(eye(ndim) + xp(i,j)*A + yp(i,j)*B)) >= sum(abs(eye(ndim) + xp(i,j)*A + yp(i,j)*B),2)) && ...
                all(diag(eye(ndim) + xp(i,j)*A + yp(i,j)*B) >= 0)
            dplot = [dplot; xp(i,j) yp(i,j)];
        end
        
    end
end

% plot feasible set with psd 
[k,~] = convhull(pplot);
plot(pplot(k,1),pplot(k,2))
hold on
% plot feasible set with DD cone
[k,~] = convhull(dplot);
plot(dplot(k,1),dplot(k,2))

