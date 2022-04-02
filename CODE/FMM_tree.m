function [ielem,itree,loct,numt,ifath,level,lowlev] = FMM_tree(x,N,maxl,levmx,ncellmx,nleafmx)
% FMM_tree: Function that creates the quad-tree structure for the Fast
%           Multipole Method solution 

% INPUTS
% --- discretization:
% - x       : coordinates of the nodes 
% - N       : number of nodes 
% --- parameters:
% - maxl    : max number of elements in a leaf 
% - levmx   : max number of level
% - ncellmx : max number of cells
% - nleafmx : max number of leaves 

% OUTPUTS 
% - ielem   : from number in tree to original number of elem 
%           (nodes ordered wrt the tree and not wrt the boundary) 
% - itree   : from global cell number to local cell number (to 0,1,2,3,...)
%           (cell number in each level with row numbering from lower left)
% - loct    : index of the 1st element in the i-th cell in the ielem vector 
%           (location of the cells in the ielem nodes of the tree vector')
% - numt    : number of elements in each cell 
% - ifath   : vector of the parent cell number 
% - level   : first cell number in a given level 

%% Parameters 

xmin = min(x,[],2);  % smallest coordinates (x and y, 2x1)
xmax = max(x,[],2);  % biggest  coordinates (x and y, 2x1)

ceilfix = @(x)ceil(abs(x)).*sign(x);

xmin = ceilfix(xmin);
xmax = ceilfix(xmax);

%% Zero Level 

% Nodes are initially ordered 
ielem = (1:N)';

% Initialization 
itree(1) = 0;   % the first cell is of class zero in its level 
level(1) = 1;   % the first cell of level zero is the cell number one 
level(2) = 2;   % the first cell of level one is the cell number two 
loct(1)  = 1;   % the first node of cell one is the first of ielem 
ifath(1) = 1;   % the parent cell of cell one is one itself
numt(1)  = N;   % cell one contains all the nodes 

ndivx    = 1;   % number of divisions (zero in this case, *2 everytime)
lowlev   = 1;   % lowest level in the tree 
nleaf    = 0;   % number of leaves 
% nswa     = 0;   % number of 
% iwork(1) = 0;   % we'll later store the number of leaves in this variable 

%% All other levels 

lev = 1;

while lev < levmx    % for every new level 
        
    % If all leaves at previous iteration end cicle
    if level(lev+1)==level(lev)   %no leaves were added
        break           
    end
    
    % Next level starts from the first available cell number 
    level(lev+1+1) = level(lev+1); 
     
    ndivxp = ndivx;      % number of divisions up to previous level
    ndivx  = 2*ndivxp;   % number of divisions up to current level (+1)
    
    dxp = (xmax-xmin)/ndivxp;    % current cell size (in x and in y)

    for cell = level(lev):level(lev+1)-1     % for every cell of the prev level
        
        % local cell number in the current level
        cell_l = itree(cell);       
        
        % if it was not a leaf 
        % if the number of nodes in the cell is still greater than the threshold
        if numt(cell) > maxl || (lev <= 2 && numt(cell)~=0) 
            
            % x coord of cell_l in its level
            itrpx = mod(cell_l,ndivxp);  
            % y coord of cell_l in its level
            itrpy = floor(cell_l/ndivxp);        
            % coordinates for separation
            xsep  = xmin + ([itrpx;itrpy] + 0.5).*dxp;       
            
            % Separates the nodes in the cell wrt y
            [ielem,nsepy ] = bisec(x,loct(cell)        ,numt(cell)      ,xsep(2),2,ielem);
            % Separates the nodes in the cell wrt x (lower part)
            [ielem,nsepx1] = bisec(x,loct(cell)        ,nsepy           ,xsep(1),1,ielem);
            % Separates the nodes in the cell wrt x (upper part)
            [ielem,nsepx2] = bisec(x,loct(cell)+nsepy  ,numt(cell)-nsepy,xsep(1),1,ielem);
            
            % NB. ielem has now been reordered in the part concerning cell
            %     we also have informations on how many nodes in each child cell
            nwk = zeros(4,1);
            % Number of nodes of each child 
            nwk(1) = nsepx1;                            % number of nodes in the child 0
            nwk(2) = nsepy - nsepx1;                    % number of nodes in the child 1
            nwk(3) = nsepx2;                            % number of nodes in the child 2
            nwk(4) = numt(cell) - nsepy - nsepx2;       % number of nodes in the child 3 
            
            locc = loct(cell);   % index of the 1st elem of the current cell in the elem vector 
            
            % Numbering and adding the new children cells to the quad tree
            for icld_y = 0:1            % child index wrt y    
                for icld_x = 0:1          % child index wrt x
                    icld = icld_y * 2 + icld_x + 1;  % child local index (+1)
                    % If child has at least one node 
                    if nwk(icld) > 0         
                        % Set nrel as the first avalable cell number 
                        nrel = level(lev+1+1);    
                    
                        % Warning 
                        if nrel > ncellmx              % if nmax of cells is reached 
                            disp(" ncellmx Error: reached max number of cells")
                        end

                        % Add new cell to the tree
                        itree(nrel) = ((itrpy*2+icld_y)*ndivxp + itrpx)*2 + icld_x;  % local cell numbering
                        % Store its first elem 
                        loct(nrel)  = locc;                                          % location on ielem of nodes of nrel cell
                        % Set its number of elements 
                        numt(nrel)  = nwk(icld);                                     % number of elements in the new cell
                        % Save its father
                        ifath(nrel) = cell;                                          % parent cell number 
                        % Possibly update the lowest level 
                        lowlev      = max(lev,lowlev);                               % lowest level update is current level 

                        % if number of elements is lower then threshold we have a leaf!
                        if (lev~=1 && (numt(nrel) <= maxl || lev == levmx))
                            nleaf = nleaf + 1; 
                            
                            if nleaf > nleafmx
                                disp(" nleafmx error: max number of leaves reached")
                            end
                        end
                        % Change the next available cell number 
                        level(lev+1+1) = nrel + 1;
                        % Location of current child elemts is adjourned 
                        locc = locc + nwk(icld);    
                    end
                end
            end 
            
        end 
    end 
    lev = lev +1;
end


end
