function [ax] = matvec(N,u,lowlev,x,y,node,dnorm,bc,...
                       ielem,itree,level,loct,numt,ifath,...
                       nexp,ntylr,maxl)
% Matvec : function that computes the matrix vector multiplication
%          exploiting the Fast Multipole Method 

% INPUTS 
% ------------------------------------------------- regarding info and data 
% - N      : number of nodes 
% - u      : vector of nodal values of the solution at last step 
% - lowlev : lowest level in the tree
% - x      : nodes
% - y      : endpoints coordinates 
% - node   : element connectivity
% - dnorm  : normal vector to every element
% - bc     : boundary condition type in each node
% --------------------------------------- regarding the quad-tree structure
% - ielem  : from number in tree to original number of elem 
% - itree  : from global cell number to local cell number (to 0,1,2,3,...)
% - loct   : index of the 1st element in the i-th cell in the ielem vector 
% - numt   : number of elements in each cell 
% - ifath  : parent cells numbers 
% - level  : first cell number in a given level 
% --------------------------------------- 
% - nexp   : truncation order of the moment expansion 
% - ntylr  : truncation order of the local expansion 
% - maxl   : max number of elements in 

%
% OUTPUT
% - ax     : vector of the values of the matrix vector multiplication  

%% Parameters 

xmin = min(x,[],2);  % smallest coordinates (x and y, 2x1)
xmax = max(x,[],2);  % biggest  coordinates (x and y, 2x1)
ceilfix = @(x)ceil(abs(x)).*sign(x);
xmin = ceilfix(xmin);
xmax = ceilfix(xmax);

%% Downward pass : compute the moments 
% Compute the moments of the multipole expansion 

[a] = upward(u,y,node,dnorm,bc,xmax,xmin,...
             ielem,itree,level,loct,numt,ifath,...
             nexp,lowlev,maxl);


%% Upward pass : compute the local coefficients 
% Exploits the moment of the multiple expansions computed before to obtain
% the coefficients of the local expansion and the matrix vector
% multiplication Ax with the solution at the previous level

[ax,~] = downward(u,N,x,y,node,dnorm,bc,a,xmax,xmin,...
                  ielem,itree,level,loct,numt,ifath,...
                  nexp,ntylr,lowlev,maxl);


end