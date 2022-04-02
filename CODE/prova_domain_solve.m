clc
clear
close all

% Boundary Discretization and Data
problem_prova2

% Exact Solution 
u_ex = [377.258*ones(N,1);-400*ones(N,1)];

% My solution 
u    = [277.258*ones(N,1);-400*ones(N,1)];


% Domain Discretization 
w = 10;
R = linspace(1,2,w);
R(1)   = R(1)+0.01;
R(end) = R(end)-0.01;

% Solution inside circular domain 
xfield = [];
for i = 1:w
    xfield = [xfield ,[R(i)*cos(dalpha*(0:N-1));R(i)*sin(dalpha*(0:N-1))]];
end
nfield = size(xfield,2);

f = domain_solve(nfield,xfield,2*N,y,bc,node,dnorm,u);

sol = f(1:N:end);

% Display solution 
disp('Sol values are :')
disp(sol')

% Plot
plot(xfield(1,:),xfield(2,:),'c--')
plot(xfield(1,:),xfield(2,:),'dc')
legend('Boundary Nodes','Int Boundary','Ext Boundary','Discretiz Levels','Discretiz Nodes')

