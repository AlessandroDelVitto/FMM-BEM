function f = domain_solve(nfield,xfield,n,y,bc,node,dnorm,u)
% Function that solves the BIE in the field points xfield inside the domain
% given the boundary solution 

% INPUTS :
%  - nfield : number of field points
%  - xfield : coordinateds of the field points 
%  - n      : number of nodes 
%  - y      : element endpoints coordinates
%  - bc     : boundary conditions type and value
%  - node   : connectivity matrix 
%  - dnorm  : normal versor in the boundary nodes 
%  - u      : solution at nodes on boudary 
% OUTPUT :
%  - f      : solution in field points 

pi2 = pi*2;

% Initialization 
f = zeros(nfield,1);

% Loop over all elements//nodes 
for j = 1:n 
    
    % If Dirichlet BC 
    if bc(1,j)==1
        f0 = bc(2,j);
        df0 = u(j);
    % If Neumann BC 
    elseif bc(1,j)==2
        f0 = u(j);
        df0 = bc(2,j);
    end
    
    % Element length
    len = sqrt((y(1,node(2,j)) - y(1,node(1,j)))^2 + ...
         (y(2,node(2,j)) - y(2,node(1,j)))^2);
     
    % Loop over all field points inside the domain
    for i = 1:nfield 
        % Distance vector wrt element j endpoints
        x11 = y(1,node(1,j)) - xfield(1,i);
        x21 = y(2,node(1,j)) - xfield(2,i);
        x12 = y(1,node(2,j)) - xfield(1,i);
        x22 = y(2,node(2,j)) - xfield(2,i);
        % Distance abs value wrt element j endpoints
        r1 = sqrt(x11^2 + x21^2);
        r2 = sqrt(x12^2 + x22^2);
        % Distance wrt tangent to the vector 
        d = x11*dnorm(1,j) + x21*dnorm(2,j);
        % Distance of endpoints wrt reference level
        t1 = -x11*dnorm(2,j) + x21*dnorm(1,j);
        t2 = -x12*dnorm(2,j) + x22*dnorm(1,j);
        % Angle difference
        ds = abs(d);
        theta1 = atan2(t1,ds);
        theta2 = atan2(t2,ds);
        dtheta = theta2 - theta1;
        
        aa = (-dtheta*ds + len + t1*log(r1) - t2*log(r2))/pi2;
        
        if d < 0 
            dtheta = -dtheta;
        end
        
        bb = -dtheta/pi2;
        
        % INFLUENCE OF NODE J ON THE FIELD POINT I 
        f(i) = f(i) + aa*df0 - bb*f0;
        
    end
end

end


