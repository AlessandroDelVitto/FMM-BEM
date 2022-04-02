function b = fmm_b(N,x,y,node,dnorm,bc,xmax,xmin,...
                   ielem,itree,level,loct,numt,ifath,...
                   nexp,ntylr,lowlev,maxl)

% Switch the BC type
bc_new = bc;
bc_new(1,bc(1,:)==1)=2;
bc_new(1,bc(1,:)==2)=1;
bc = bc_new;

% The BC are ordered accordingly to the tree in the sol u
u  = bc(2,ielem);

[a] = upward(u,y,node,dnorm,bc,xmax,xmin,...
             ielem,itree,level,loct,numt,ifath, ...
             nexp,lowlev,maxl);

[ax,~] = downward(u,N,x,y,node,dnorm,bc,a,xmax,xmin,...
                  ielem,itree,level,loct,numt,ifath,...
                  nexp,ntylr,lowlev,maxl);
 
% Return 
b = -ax;

end
