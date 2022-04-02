% To verify the code functioning I take the simple problem in the domain
% omega = [0,1]^2 in R^2 with homogeneous Dirichlet boundary conditions 

% BEM DISCRETIZATION
% I have to define 
% x,y,N,node,dnorm 

% element endpoints
len = 1;
resolution = 20;
dx = len/resolution;
dy = len/resolution;
endpoints_x = 0:dx:len;
endpoints_y = 0:dy:len;
n = numel(endpoints_x);

y = [[endpoints_x;zeros(1,n)],[ones(1,n-1);endpoints_y(2:end)],[flip(endpoints_x(1:end-1));ones(1,n-1)],[zeros(1,n-2);flip(endpoints_y(2:end-1))]];

% nodes 
nodes_x = dx/2:dx:1-dx/2;
nodes_y = dx/2:dx:1-dx/2;
nodes_n = numel(nodes_x);
x = [[nodes_x;zeros(1,nodes_n)],[ones(1,nodes_n);nodes_y(1:end)],[flip(nodes_x(1:end));ones(1,nodes_n)],[zeros(1,nodes_n);flip(nodes_y(1:end))]];

% number of elements 
N = length(x);

% element connectivity 
node = zeros(2,N);
node(1,:)   = 1:N;
node(2,:)   = node(1,:)+1;
node(2,end) = 1;

% dnorm 
dnorm = zeros(2,N);
dnorm(2,1:N/4)       = -1;
dnorm(1,N/4+1:N/2)   =  1;
dnorm(2,N/2+1:N/4*3) =  1;
dnorm(1,N/4*3+1:N)   = -1;
for i = 1:N
    h1 = y(2,node(2,i)) - y(2,node(1,i));
    h2 = -y(1,node(2,i)) + y(1,node(1,i));
    el = sqrt(h1^2 + h2^2);
    dnorm(1,i) = h1/el;
    dnorm(2,i) = h2/el;
end

% boundary condition (Dirichlet Type == 1)
bc = [ones(1,N);ones(1,N)];

% plot 
figure()
plot(x(1,:),x(2,:),'dr')
hold on 
plot(y(1,:),y(2,:),'ob')
