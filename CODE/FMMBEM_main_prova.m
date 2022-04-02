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

% 1) Validation problem #1
problem_prova


%% TREE 

% Number of elements 
N = size(x,2);

% Parameters 
maxl    = 5;     % max number of elements in a leaf
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

% 1) It is completely null in the Homogeneous-Dirichlet case 
b = zeros(N,1);

% 2) For the second problem see paper 
xmax = max(x,[],2);
xmin = min(x,[],2);
ceilfix = @(x)ceil(abs(x)).*sign(x);
xmin = ceilfix(xmin);
xmax = ceilfix(xmax);
% 
% b = fmm_b(N,x,y,node,dnorm,bc,xmax,xmin,    ...
%           ielem,itree,level,loct,numt,ifath,...
%           nexp,ntylr,lowlev,maxl);
b_2 = bvector(x,y,bc,node,dnorm,N);

             
%% LINEAR SYSTEM SOLUTION 

% The function matvec exploits the FMM to compute the product A*lambda
% from the previous iteratve solution u_(k-1)
fmm_solver = @(u) matvec(N,u,lowlev,x,y,node,dnorm,bc,...
                         ielem,itree,level,loct,numt,ifath,...
                         nexp,ntylr,maxl);

% Gmres solution of the system 
restart = [];
tol     = 1e-8;
maxit   = N;
u0      = zeros(N,1);

sol = gmres(fmm_solver,b,restart,tol,maxit,[],[],u0);

%        ---> fmm_solver  handle function that computes the product A*lambda with FMM
%             b           RHS of the system 
%             sol         solution of the system ---> flux nodal values 


% Final Solution Plot 
figure()
plot(1:N,sol,'rd')
hold on 
ex_sol = ones(N,1);
plot(1:N,ex_sol(ielem),'b-')
title('Boundary Nodal Values of the FMM solution')
grid on 
xlabel('Node')
ylabel('Nodal value')
legend('FMM-BEM sol','Exact sol','Location','east')



%% WHOLE DOMAIN SOLUTION

% Solve the BIE with the known values of phi and q on the boundary


% Plot the solution 
% figure()




