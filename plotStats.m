function []  = plotStats(iter)


iterations = length(iter);

% extract data
duals           = zeros(iterations,1);
primal          = zeros(iterations,1);
cost            = zeros(iterations,1);
nabla_Langrange = zeros(iterations,1);
conVio = zeros(iterations,1);
for k = 1:iterations
    duals(k) = iter{k}.dual;
    primal(k) = iter{k}.primal;
    cost(k) = iter{k}.cost;
    nabla_Langrange(k) = iter{k}.nabla_Langrange;
    conVio(k)      = iter{k}.conVio;
end


figure(2)
plot(1:iterations,primal ,'-*')
xlabel('Iterations')
ylabel('Primal')

figure(3)
plot(1:iterations,duals,'-*')
xlabel('Iterations')
ylabel('Dual')

figure(4)
plot(1:iterations,cost,'-*')
xlabel('Iterations')
ylabel('Cost')

figure(5)
plot(1:iterations,nabla_Langrange,'-*')
xlabel('Iterations')
ylabel('\nabla L')

figure(6)
plot(1:iterations,conVio,'-*')
xlabel('Iterations')
ylabel('||\theta(x_k)||_2')

end