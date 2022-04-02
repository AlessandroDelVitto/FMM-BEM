% FOR ALL THE LEAVES 

% If icell is a leaf or if we're in the last level
if numt(icell) <= maxl || lev == lowlev
   fact = 1;
   % Multiply the local coeff by 1/order
   for itylr = 1:ntylr
        fact = fact/itylr;
        b(itylr+1,icell) = b(itylr+1,icell)*fact;
   end
   % For every node in the cell 
   for in = 1:numt(icell)
       % Element number in the ielem vector 
       inax = loct(icell)+in-1; 
       % Original element number
       indx = ielem(inax); 
       % Local coeff of max order 
       zp = b(ntylr+1,icell);
       % Vector from the center of icell to the node inside the cell
       z0 = complex(x(1,indx)-c(1),x(2,indx)-c(2));
       % For every order starting from the last one -1 backwds
       for itylr = ntylr:-1:1
           % Local expansion
            zp = zp*z0 + b(itylr,icell);  
       end
       zp = zp/(pi*2);
       % Array Ax
       ax(inax) = ax(inax) + real(zp); 
   end
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% FOR ALL THE LEAVES 

% If icell is a leaf or if we're in the last level
if numt(icell) <= maxl || lev == lowlev
   % For every node in the cell 
   for in = 1:numt(icell)
       % Element number in the ielem vector 
       inax = loct(icell)+in-1; 
       % Original element number
       indx = ielem(inax); 
       % Vector from the center of icell to the node inside the cell
       z0 = complex(x(1,indx)-c(1),x(2,indx)-c(2));
       
       zp = 0;
       for l = 0:ntylr
           % Local expansion
            zp = zp + 1/(2*pi)*b(l+1,icell)*(z0)^l/factorial(l);
       end
       
       % Array Ax
       ax(inax) = ax(inax) + real(zp); 
   end
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% FOR ALL THE LEAVES 

% If icell is a leaf or if we're in the last level
if numt(icell) <= maxl || lev == lowlev
   % For every node in the cell 
   for in = 1:numt(icell)
       % Element number in the ielem vector 
       inax = loct(icell)+in-1; 
       % Original element number
       indx = ielem(inax); 
       % Vector from the center of icell to the node inside the cell
       z0 = complex(x(1,indx)-c(1),x(2,indx)-c(2));
       
       I = [1;cumprod(ones(ntylr,1)*z0)];
       fact = 2*pi*[1 cumprod(1:ntylr)]';
       I = I ./ fact ;
       zp = b(:,icell) * I ;
       
       % Array Ax
       ax(inax) = ax(inax) + real(zp); 
   end
end
