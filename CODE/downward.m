function [ax,b] = downward(u,N,x,y,node,dnorm,bc,a,xmax,xmin,...
                           ielem,itree,level,loct,numt,ifath,...
                           nexp,ntylr,lowlev,maxl)
% downward : function that computes the coefficients of the local expansions
%            of all the nodes in the tree from the level two cells to the
%            leaves of the quad-tree structure. 
%            1) Cells in the interaction list of C are dealt with M2L translation (eq 25)
%            2) far cells form C are dealt with L2L translations (eq 27 & 28).
%
% INPUTS 
% ------------------------------------------------- regarding info and data 
% - u      : solution at previous step 
% - N      : number of nodes 
% - y      : coordinates of the endpoints of each element 
% - node   : element connectivity 
% - dnorm  : normal vector to every element
% - bc     : boundary condition type in each node
% - a      : moment matrix 
% - xmax   : greatest coordinates
% - xmin   : smallest cordinates 
% --------------------------------------- regarding the quad-tree structure
% - ielem  : from number in tree to original number of elem 
%            nodes ordered accordingly to the tree 
% - itree  : from global cell number to local cell number (to 0,1,2,3,...)
%            cell numbers in each level
% - loct   : index of the 1st element in the i-th cell in the ielem vector 
%            location of the cells in the ielem nodes of the tree vector 
% - numt   : number of elements in each cell 
% - ifath  : parent cells numbers 
% - level  : first cell number in a given level 
% ------------------------------------- regarding the parameters of the FMM
% - nexp  : number of coefficients in the moment expansion
% - ntylr : number of terms in local expansion 
% % - rwork : vector of block preconditioning matrix values 
% % - iwork : vector informations about the block preconditioning matrix 
% %           location and dimension of each single block
%
% OUTPUTS 
% - ax    : resulting vector of multiplications Ax
% - b     : local expansion coefficients of al the cells (matrix)
%
% NB. a cell C of level two has no far cells so it just has M2L
%     translations for the cells in the interaction listof C

% Inizialization 
b = zeros(ntylr+1,level(lowlev+1+1)-1);

% Vector of the results of the multiplication Ax (n terms)
ax = zeros(N,1);

% Leaves count 
leaf_n = 0;

% For every level starting from level 2 (downward movement)
for lev = 2:lowlev   % 110
    % Division factor from the original cell to the cells in lev
    ndivx = 2^lev;
    % Size of a cell of the current level lev 
    dx = (xmax-xmin)/ndivx;
    % Loop for all cells in level lev 
    for icell = level(lev+1):level(lev+1+1)-1  % 120
        % Local index of the current cell {0,1,2,...}
        itr = itree(icell);
        % Position of the cell in the current level lev
        itrx = mod(itr,ndivx);
        itry = floor(itr/ndivx); 
        % Center of the cell
        c = xmin + ([itrx;itry] + 0.5) .* dx;
        % Parent cell position in its level (lev-1)
        itrxp = floor(itrx/2);
        itryp = floor(itry/2);  
        
        %% CONTRIBUTION OF FAR CELLS FROM C (through parent cells)
        
        % If not level 2 you have far cells to deal with 
        if lev ~= 2 
            % Center of the parent cell of icell
            cp = xmin + ([itrxp;itryp]+0.5) .* (dx*2);
            % Vector from center of parent to center of icell
            z0 = complex(c(1)-cp(1), c(2)-cp(2));
            % Global cell # of the parent cell
            io = ifath(icell); 
            
            % L2L TRANSLATION from parent cell
            % NB the l-th local term in the expansion of the child cell is
            %    influenced only by higher order local terms of the parent
            %    cell expansion
            
            % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MINE 
            % Compute the influence of the ntylr coeff of the parent cell
            % on the lower degree coeff of the child cell
%             for l = 0:ntylr
%                 zi = 1;
%                 % Each influence only the lower degrees
%                 for m = l:ntylr
%                     % Sum to the icell local exp coefficient of order m
%                     % the influence of the local coeff of order k+m of the parent cell
%                     b(l+1,icell) = b(l+1,icell) + zi * b(m+1,io);
%                     zi = zi*z0/(l+1);
%                 end
%             end
            
            % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> HIS
            zi = complex(1,0);
            % Compute the influence of the ntylr coeff of the parent cell
            % on the lower degree coeff of the child cell
            for k = 0:ntylr
                % Each influence only the lower degrees
                for m = 0:ntylr-k
                    % Sum to the icell local exp coefficient of order m
                    % the influence of the local coeff of order k+m of the parent cell
                    b(m+1,icell) = b(m+1,icell) + zi * b(k+m+1,io);
                end
                zi = zi*z0/(k+1);
            end
            % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            
        end
        
        %% CONTRIBUTION OF CELLS IN THE INTERACTION LIST OF ICELL OR ADJACENT
        
        % For every other cell of the level lev
        for jcell = level(lev+1):level(lev+1+1)-1  %130
            % Local index of the current cell
            jtr  = itree(jcell);
            % Local position of the current cell
            jtrx = mod(jtr,ndivx);
            jtry = floor(jtr/ndivx);
            % Local position of the parent cell
            jtrxp = floor(jtrx/2);
            jtryp = floor(jtry/2);
            
            % IF JCELL IS NOT FAR CELL FROM ICELL
            
            % If parents are neighbours then jcell is in the interaction list 
            if abs(itrxp-jtrxp) <= 1 && abs(itryp-jtryp) <= 1
                
                % a) If cells are NOT neighbours use M2L translation
                if abs(itrx-jtrx)>1 || abs(itry-jtry)>1
                    % Center of the j cell
                    cc = xmin + ([jtrx;jtry] + 0.5) .* dx;
                    % Vector from the center of the interaction list cell jcell to the cell icell 
                    z0 = complex(c(1)-cc(1), c(2)-cc(2));
                   
                    % >>>>>>>>>>>>>>>>>>>>>>>>>>> MINE
%                     % Local coeff of order zero from moment of order zero 
%                     b(1,icell) = b(1,icell) - log(z0) * a(1,jcell);
%                     
%                     for l = 0:ntylr
%                         sgn = (-1)^l;
%                         % k=0 term
%                         if l>0
%                             O_lk = factorial(l-1)/z0^(l);
%                             b(l+1,icell) = b(l+1,icell) + sgn*O_lk*a(1,jcell); 
%                         end
%                         % k>0 terms
%                         for k = 1:nexp
%                             % M2L formula
%                             O_lk = factorial(k+l-1)/z0^(l+k);
%                             b(l+1,icell) = b(l+1,icell) + sgn*O_lk*a(k+1,jcell); 
%                         end
%                     end
                    
                    % >>>>>>>>>>>>>>>>>>>>>>>>>>> HIS
                    % Local coeff of order zero from moment of order zero 

                    b(1,icell) = b(1,icell) - log(z0) * a(1,jcell);

                    % First term
                    zo = complex(1,0);

                    % Compute the influence of every moment of jcell on the
                    % proper local coeff of icell
                    for m = 1:nexp+ntylr
                        % Divide by the distance between the centers 
                        zo = zo/z0;
                        % Compute the coeff span that is influenced by Mm
                        kmin = max(0,m-nexp);
                        kmax = min(m,ntylr);
                        % Compute sign of kmin term
                        sgn = (-1)^kmin;
                        % For every term from kmin to kmax 
                        for k = kmin:kmax
                            % M2L formula 
                            b(k+1,icell) = b(k+1,icell) + sgn*zo*a(m-k+1,jcell); 
                            % Change sign 
                            sgn = -sgn;
                        end
                        % Multiply for the factorial factor
                        zo = zo*m;
                    end
                    % >>>>>>>>>>>>>>>>>>>>>>>>>>>
                    
                % b) If jcell and icell are adjacent then compute the direct influence (BEM) 
                % If jcell or icell are adjacent leaves 
                elseif numt(jcell)<=maxl || numt(icell)<=maxl || lev==lowlev
                    
                    % If they're re same cell
                    if icell==jcell
                        % Count that you passed in this leaf 
                        leaf_n = leaf_n+1;
                    end

                    % Compute the direct influence of jcell to icell
                    % Direct computation of the product Ax
                    ax(loct(icell):loct(icell)+numt(icell)-1) = ax(loct(icell):loct(icell)+numt(icell)-1) + ...
                          direct(ielem,loct,numt,icell,jcell,node,x,y,...
                                 dnorm,u,bc);
                end
            end
                
        end % line 130 fortrun 


        % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        % FOR ALL THE LEAVES 
        % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MINE
        % ---> compute the matrix vector multiplication 

%             % If icell is a leaf or if we're in the last level
%             if numt(icell) <= maxl || lev == lowlev
%                % For every node in the cell 
%                for in = 1:numt(icell)
%                    % Element number in the ielem vector 
%                    inax = loct(icell)+in-1; 
%                    % Original element number
%                    indx = ielem(inax); 
%                    % Vector from the center of icell to the node inside the cell
%                    z0 = complex(x(1,indx)-c(1),x(2,indx)-c(2));
% 
%                    I = [1;cumprod(ones(ntylr,1)*z0)];
%                    fact = 2*pi*[1 cumprod(1:ntylr)]';
%                    I = I ./ fact ;
%                    zp = b(:,icell)' * I ;
% 
%                    % Array Ax
%                    ax(inax) = ax(inax) + real(zp); 
%                end
%             end

        % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> HIS
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
               dz = complex(x(1,indx)-c(1),x(2,indx)-c(2));
               % For every order starting from the last one -1 backwds
               for itylr = ntylr-1:-1:0
                   % Local expansion
                    zp = zp*dz + b(itylr+1,icell);  
               end
               zp = zp/(pi*2);
               % Array Ax
               ax(inax) = ax(inax) + real(zp); 
           end
        end

        % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    end  %120
end  % 110

end