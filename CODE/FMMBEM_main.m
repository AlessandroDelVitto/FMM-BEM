clear 
close all
clc

% -------------------------------------------------------------------------
%   Alessandro Del Vitto 
%   Politecnico di Milano 
%   Mathematical Engineering 2021-2022
%   Low Frequency Computational Electromagnetics 
%   " Fast Multipole Method implementation of a 2D Boundary Element Method 
%     Code for the solution of the Laplace equation "
% -------------------------------------------------------------------------


%% DATA 
% Read data of the problem 
% Read parameters for the FMM & BEM 

% must obtain x,y,N,node,dnorm;

% 2) Validation problem #2
problem_prova2


%% TREE 

% Number of elements 
N = size(x,2);

% Parameters 
maxl    = 20;     % max number of elements in a leaf
levmx   = 10;    % max number of levels in the tree 
ncellmx = 100;   % max number of cells in the tree 
nleafmx = 50;    % max number of leaves in the tree

% Construct the quad-tree structure 
[ielem,itree,loct,numt,ifath,level,lowlev] = FMM_tree(x,N,maxl,levmx,ncellmx,nleafmx);


%% RIGHT HAND SIDE 

% PARAMETERS of the FMM
% Truncation order of the moment expansion
nexp  = 15;
% Truncation order of the local expansion
ntylr = 15;

% 2) For the second problem (see paper)
xmax = max(x,[],2);
xmin = min(x,[],2);
ceilfix = @(x)ceil(abs(x)).*sign(x);
xmin = ceilfix(xmin);
xmax = ceilfix(xmax);

b = fmm_b(N,x,y,node,dnorm,bc,xmax,xmin,    ...
          ielem,itree,level,loct,numt,ifath,...
          nexp,ntylr,lowlev,maxl);
b_2 = bvector(x,y,bc,node,dnorm,N);

% Plot RHS difference 
figure()
plot(b,'r-o')
hold on 
grid on 
plot(b_2(ielem),'b-o')
title('RHS FMM approximation')
legend('FMM-RHS','BEM-RHS','Location','East')

figure()
plot(b-b_2,'r-o','Linewidth',1)
grid on 
title('RHS approximation error with FMM')

             
%% LINEAR SYSTEM SOLUTION 

% The function matvec exploits the FMM to compute the product A*lambda
% from the previous iteratve solution u_(k-1)
fmm_solver = @(u) matvec(N,u,lowlev,x,y,node,dnorm,bc,...
                         ielem,itree,level,loct,numt,ifath,...
                         nexp,ntylr,maxl);

% Gmres solution of the system 
restart = [];
tol     = 1e-8;
maxit   = 72;
u0      = [300*ones(N/2,1);-400*ones(N/2,1)];

sol = gmres(fmm_solver,b_2,restart,tol,maxit,[],[],u0);

%        ---> fmm_solver  handle function that computes the product A*lambda with FMM
%             b           RHS of the system 
%             sol         solution of the system ---> flux nodal values 


% Final Solution Plot 
figure()
plot(1:N,sol,'rd')
hold on 
ex_sol = [377.258*ones(N/2,1);-400*ones(N/2,1)];
plot(1:N,ex_sol(ielem),'b-')
title('Boundary Nodal Values of the FMM solution')
grid on 
xlabel('Node')
ylabel('Nodal value')
legend('FMM-BEM sol','Exact sol','Location','east')



%% WHOLE DOMAIN SOLUTION

% Solve the BIE with the known values of phi and q on the boundary

% Reorder solution 
sol(ielem) = sol;

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

f = domain_solve(nfield,xfield,N,y,bc,node,dnorm,sol);

int_sol = f(1:N:end);

% Display solution 
disp('Sol values are :')
disp(int_sol')

% Plot
figure()
plot(x(1,:),x(2,:),'dr')
hold on 
plot(y1(1,:),y1(2,:),'ob-')
plot(y2(1,:),y2(2,:),'ob-')
grid on 
xlabel('x')
ylabel('y')
title('Discretized domain \omega')
plot(xfield(1,:),xfield(2,:),'c--')
plot(xfield(1,:),xfield(2,:),'dc')
legend('Boundary Nodes','Int Boundary','Ext Boundary','Discretiz Levels','Discretiz Nodes')




