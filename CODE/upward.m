function [a] = upward(u,y,node,dnorm,bc,xmax,xmin,     ...
                      ielem,itree,level,loct,numt,ifath, ...
                      nexp,lowlev,maxl)
% Upward : function that computes the moments on all levels starting from 
%          the leaves up to the second level. The moments are approximated
%          up to p terms. Leaves have a direct computation whereas parent
%          cells' moments are computed by an M2M translation from children.
% 
% INPUTS 
% ------------------------------------------------- regarding info and data 
% - u      : solution at previous step 
% - y      : coordinates of the endpoints of each element 
% - node   : element connectivity 
% - dnorm  : normal vector to every element
% - bc     : boundary condition type in each node
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
% -------------------------------------------------- parameters of the FMM
% - nexp   : number of terms in the multipole expansion
% - lowlev : lowest level in the tree 
% - maxl   : max number of elements in a leaf 

% OUTPUTS
% - a      : multipole expansion moments of every cell (matrix nexp x ncells)

% Initialization 
ncells = level(lowlev+1+1)-1;       % total number of cells 
a      = zeros(nexp+1,ncells);

% Loop from lowest level to level 2 cells (Upward movement)
for lev = lowlev:-1:2           
    % Dividing factor of the considered level lev
    ndivx = 2^lev;    %--> number of cells on each row/col
    % Determine cell size in this level lev
    dx = (xmax-xmin)/ndivx;        
    % For all the cells of level lev 
    for icell = level(lev+1):level(lev+1+1)-1   
        % Local cell number (0,1,...)
        local_cell_index = itree(icell);
        % Local cell coordinates (Position of the cell in the level)
        local_cell_coord_x = mod(local_cell_index,ndivx);
        local_cell_coord_y = floor(local_cell_index/ndivx);   
        % Center of the cell
        c = xmin + ([local_cell_coord_x;local_cell_coord_y] + 0.5) .* dx;
        
        %% If leaf -> use eq 20
        if(numt(icell) <= maxl || lev == lowlev)  
            % Number of elements in the leaf 
            num = numt(icell);
            % Index of the first elem of the leaf in the ielem vector
            first = loct(icell);
            
            % Compute the k moments of the leaf 
            a(:,icell) = FMM_moment(a(:,icell),y,node,ielem,num,nexp,c,u,bc,dnorm,first);

        end
        
        %% If not leaf -> use parent cells
        if lev~=2
            % Do M2M translation to form moments of the parent cell
            
            % Center of parent cell
            cp = xmin + (floor([local_cell_coord_x;local_cell_coord_y]/2) + 0.5) .* dx * 2;
            
            % Complex vector (distance between centers) 
            % from parent cell to current cell icell
            dz = complex(c(1)-cp(1), c(2)-cp(2));
            
            % Global cell number of parent cell
            io = ifath(icell);
            
            I_kl = complex(1,0);
            % NB the m-th moment of the parent cell is influenced only by
            %    the first m moments of the children cells
            
            % For every term of the child expansion
            for k = 0:nexp 
                % Compute its influence on the higher order moments of the
                % parent cell expansion 
                for m = k:nexp
                    % M2M translation 
                    % Sum to the moment of order m of the parent cell the influence 
                    % of the moments of order m-k of the child cell icell
                    a(m+1,io) = a(m+1,io) + I_kl * a(m-k+1,icell);
                end
                % Function I_{k-l}(z_c-z'_c)
                I_kl = I_kl * dz/(k+1);
            end
            
        end
    end
end

end