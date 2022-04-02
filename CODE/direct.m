function [ax] = direct(ielem,loct,numt,icell,jcell,node,x,y,dnorm,u,bc)
% direct : Function that computes the integrals in a direct way.
%          It computes the influence of all the elements in the neighbouring
%          cell jcell to all the elements in the cell icell. (BEM EQ)
%
% INPUTS
% - inod  : vector of indexes of the nodes of icell
% - jnod  : vector of indexes of the nodes of jcell
% - node  : connectivity vector 
% - x     : nodes of the boundary discretization (coordinates vector)
% - y     : endpoints of all the elements (coordinates vector)
% - ni    : number of nodes in icell
% - nj    : number of nodes in jcell
% - dnorm : vector with the normal versor of each element 
% - indx_first_elem_j : index of the first elem of jcell in u vector 
% - u     : solution at previous step 
% % - icell : cell number of the C cell icell
% % - jcell : cell number of the neighbouring cell jcell
% % - amat  : diagonal block preconditioning matrix
% % - isw   : information about the precond matrix in the current cell
% - bc    : matrix with the BC type in each node and value on second row
%
% OUTPUTS
% - ax    : vector with the product Ax regarding icell

% Number of nodes 
ni = numt(icell);
nj = numt(jcell);

% Initialization of the influence vector
ax = zeros(ni,1);

% For every node in jcell
for j = 1:nj
    % Global index of the j-th node in jcell
    j_ind = ielem(loct(jcell)+j-1);
    % Element length
    len = sqrt((y(1,node(1,j_ind))-y(1,node(2,j_ind)))^2 +(y(2,node(1,j_ind))-y(2,node(2,j_ind)))^2);
    % For every node in icell
    for i = 1:ni
        % Global index of the i-th node in icell
        i_ind = ielem(loct(icell)+i-1);
        
        % Distance between the i-th node of icell with the end nodes of the
        % element of the j-th node of jcell
        x11 = y(1,node(1,j_ind)) - x(1,i_ind);
        x21 = y(2,node(1,j_ind)) - x(2,i_ind);
        x12 = y(1,node(2,j_ind)) - x(1,i_ind);
        x22 = y(2,node(2,j_ind)) - x(2,i_ind);
        
        % Modulus of distance with node1
        r1 = sqrt(x11^2 + x21^2);
        % Modulus of distance with node2
        r2 = sqrt(x12^2 + x22^2);
        
        % Scalar product of dist wrt node 1 and the normal derivative
        d = x11*dnorm(1,j_ind) + x21*dnorm(2,j_ind);
        
        t1 = -x11*dnorm(2,j_ind) + x21*dnorm(1,j_ind);
        t2 = -x12*dnorm(2,j_ind) + x22*dnorm(1,j_ind);
        
        ds = abs(d);
        
        % Arctangent 
        dtheta = atan2(ds*len,ds^2+t1*t2);
        
        % Influence coeff of node nj on node ni
        % analytical expression (see BEM with constant elements)
        aa = (-dtheta*ds + len + t1*log(r1) - t2*log(r2))/(pi*2);
        
        if d<0 
            dtheta = -dtheta;
        end
        
        bb = -dtheta/(pi*2);
        
        if i_ind == j_ind 
            bb = 0.5;
%             aa = len/(2*pi)*(1-log(len/2));  % mio
        end
        
        if bc(1,j_ind) == 1      
            % If Dirichlet boundary condition
            ax(i) = ax(i) - aa * u(loct(jcell)+j-1);
            
        elseif bc(1,j_ind) == 2 
            % If Neumann boundary condition 
            ax(i) = ax(i) + bb * u(loct(jcell)+j-1);
            
        else 
            disp('Error BC not valid')
        end
        
    end
end

end