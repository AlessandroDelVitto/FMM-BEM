% To verify the code functioning I take a problem with an analytically
% known solution see paper for reference 

% BEM DISCRETIZATION
% I have to define 
% x,y,N,node,dnorm 

% Number of nodes 
N = 36;
dalpha = 2*pi/N;

% Radiuses 
R1 = 1;     % a in the paper 
R2 = 2;     % b in the paper 

% element endpoints
% y1 = [R1*cos(dalpha*(0:N-1)+5/4*pi+dalpha/2);R1*sin(dalpha*(0:N-1)+5/4*pi+dalpha/2)];
% y2 = [R2*cos(dalpha*(0:N-1)+5/4*pi+dalpha/2);R2*sin(dalpha*(0:N-1)+5/4*pi+dalpha/2)];
% y  = [y2,y1];
y1 = [R1*cos(dalpha*(0:N-1));R1*sin(dalpha*(0:N-1))];
y2 = [R2*cos(dalpha*(0:N-1));R2*sin(dalpha*(0:N-1))];
y  = [y2,y1];

% nodes 
x1 = (y1(:,1:end-1)+y1(:,2:end))/2;
x1 = [x1,(y1(:,end)+y1(:,1))/2];
x2 = (y2(:,1:end-1)+y2(:,2:end))/2;
x2 = [x2,(y2(:,end)+y2(:,1))/2];
x  = [x2,x1]; 

% element connectivity 
node = zeros(2,2*N);
node(1,:)   = 1:2*N;
node(2,:)   = node(1,:)+1;
node(2,N)   = 1;
node(2,2*N) = N+1;

% dnorm 
dnorm = [cos(dalpha*(0:2*N-1)+5/4*pi+dalpha/2);sin(dalpha*(0:2*N-1)+5/4*pi+dalpha/2)];
dnorm2 = dnorm;
for i = 1:2*N
    h1 = y(2,node(2,i))  - y(2,node(1,i));
    h2 = -y(1,node(2,i)) + y(1,node(1,i));
    el = sqrt(h1^2 + h2^2);
    dnorm(1,i) = h1/el;
    dnorm(2,i) = h2/el;
end

% boundary condition (Dirichlet Type == 1, Neumann type == 2)
bc1 = [2*ones(1,N),ones(1,N)];
bc2 = [200*ones(1,N),100*ones(1,N)];
bc  = [bc1;bc2];

% plot 
figure()
plot(x(1,:),x(2,:),'dr')
hold on 
plot(y1(1,:),y1(2,:),'ob-')
plot(y2(1,:),y2(2,:),'ob-')
grid on 
xlabel('x')
ylabel('y')
title('Discretized domain \omega')
legend('Boundary Nodes','Int Boundary','Ext Boundary')

