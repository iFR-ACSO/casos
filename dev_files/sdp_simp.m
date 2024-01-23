% test block diagonalization in SOS program

% free variable
x = casos.PS('x', 6 ,1);

% given polynomial 
f = [1 + x(1)^4 + x(1)*x(2)+x(2)^4 + x(3)^2 ;
      x(2)^2 + x(1)^2*x(3)^2;
      x(1)^2 + x(1)^4+x(1)^6;
      1 + x(1)^4*x(2)^2+x(1)^2*x(2)^4];


%+ ... x(5:6)'*diag([2 1])*x(5:6) + x(6)^6 + x(5:6).^2'*diag([2 5])*x(5:6).^2;

tic
% obtain gram basis
[Z,K_cone,z] = f.grambasis([], 1);

A = [];
b = [];
for w=1:length(K_cone)

gdiff = f(w);
i=0;
Qcoeffs = [];
for j=1:length(K_cone{w})
    Q = casos.PS.sym(['Q' num2str(j) num2str(w)], 1, [K_cone{w}(j) K_cone{w}(j)]);
    gdiff = gdiff - z{w}(i+1:i+K_cone{w}(j))'*Q*z{w}(i+1:i+K_cone{w}(j));
    i = i + K_cone{w}(j);
    Qcoeffs = [Qcoeffs, Q.coeffs];
end

[Qdiff,~]= gdiff.poly2basis(Z{w}); 

a1 = Qdiff.jacobian(Qcoeffs);
A = mdiag(A, casadi.DM(a1).full());
    
eval_function = casadi.Function('eval_function', {Qcoeffs}, {Qdiff});
b1 = eval_function(zeros(length(Qcoeffs),1));
b = [ b; -b1.full()];

end

% cone size
K.s = cell2mat(K_cone);
c = zeros(sum(cell2mat(K_cone).^2),1);
time_pre = toc;

% build SDP problem and solve using sedumi
[x,y,info] = sedumi(A,b,c,K);
fprintf('\nsedumi time: %d \n',sum(info.timing))
fprintf('pretime: %d \n',time_pre)
fprintf('Total Time: %d \n', sum(info.timing) + time_pre)   

