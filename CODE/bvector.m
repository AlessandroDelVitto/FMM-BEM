function b = bvector(x,y,bc,node,dnorm,N)

pi2 = pi*2;
u = zeros(N,1);
for j = 1:N     % Loop on field points (Column)
    al = sqrt((y(1,node(2,j)) - y(1,node(1,j)))^2 +(y(2,node(2,j)) - y(2,node(1,j)))^2);
    for i = 1:N       % Loop on source points (Row)
        % Compute parameters used in the formulas for the two intergals
        x11 = y(1,node(1,j)) - x(1,i);
        x21 = y(2,node(1,j)) - x(2,i);
        x12 = y(1,node(2,j)) - x(1,i);
        x22 = y(2,node(2,j)) - x(2,i);
        r1 = sqrt(x11^2 + x21^2);
        r2 = sqrt(x12^2 + x22^2);
        d = x11*dnorm(1,j) + x21*dnorm(2,j);
        
        t1 = -x11*dnorm(2,j) + x21*dnorm(1,j);
        t2 = -x12*dnorm(2,j) + x22*dnorm(1,j);
        ds = abs(d);
        theta1 = atan2(t1,ds);
        theta2 = atan2(t2,ds);
        dtheta = theta2 - theta1;
        aa = (-dtheta*ds + al + t1*log(r1) - t2*log(r2))/pi2;
        if d < 0
            dtheta = -dtheta;
        end
        bb = -dtheta/pi2;
        if i == j 
            bb = 0.5;
        end
        if bc(1,j) == 1 
            % Potential is given
            u(i) = u(i) - bb*bc(2,j); 
        end
        if bc(1,j) == 2 
            % Derivative is given
            u(i) = u(i) + aa*bc(2,j); 
        end
    end
end

b = u;

end