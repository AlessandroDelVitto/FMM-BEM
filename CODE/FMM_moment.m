function [a] = FMM_moment(a,y,node,ielem,num,nexp,c,u,bc,dnorm,first)
% FMM_moment : function that computes the moments of a given leaf with
%              the direct computation, formula (20)
% INPUTS
% - ielem : elements ordered accordingly to the cells of the tree 
% - num   : number of elements in the leaf 
% - dnorm : normal vector at each node 
% - bc    : bc type (1 for D, 2 for N) 
% - node  : element connectivity 
% - y     : coordinates of the endpoints of the elements 
% - c     : coordinates of the center of the cell 
% - nexp  : number of expansion terms 
% - u     : solution vector at previous step 
% - first : index of the first element of ielem in the current leaf 

% OUTPUTS
% - a     : moment matrix's column related to the leaf in consideration 

% Inizialization 
a = zeros(nexp+1,1);

% For every element in the leaf
for i = 0:num-1
    % Element global number
    nelm = ielem(first+i);
    % Numbers of the nodes at the endpoints of the element 
    n1 = node(1,nelm);
    n2 = node(2,nelm);
    % Distances from the center of the cell 
    z1 = complex(y(1,n1)-c(1), y(2,n1)-c(2));
    z2 = complex(y(1,n2)-c(1), y(2,n2)-c(2));
    
    % omega bar 
    zwbar = conj(z2-z1);
    zwbar = zwbar/abs(zwbar);
    
    % Projections 
    zp1 = z1 * zwbar;
    zp2 = z2 * zwbar;
    
    % Complex normal versor in the element node 
    znorm = complex(dnorm(1,nelm),dnorm(2,nelm));
    
    if bc(1,nelm) == 1       % BC of type one (D) 
        % Assign values to phi and q
        phi = 0;
        q = u(first+i);
    elseif bc(1,nelm) == 2   % BC of type two (N)
        % Assign values to phi and q
        phi = u(first+i);
        q = 0;
    end
    
    % Moment of order zero 
    a(1) = a(1) - (zp2-zp1)*q;
    
    % Moments of order up to nexp
    for k = 1:nexp
        
        % F kernel
        a(k+1) = a(k+1) + (zp2-zp1) * znorm * phi; 

        % G kernel
        zp1 = zp1 * z1 /(k+1);
        zp2 = zp2 * z2 /(k+1);
        
        a(k+1) = a(k+1) - (zp2-zp1)*q;   
    end
end

end