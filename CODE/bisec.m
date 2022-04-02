function [ielem,nsep] = bisec(x,ielem_first,n,xsep,dim,ielem)
% bisec: Function that separates the nodes of a given cell in the specified
%        dimension wrt a given separation coordinate xsep 

% INPUTS 
% - x           : nodes coordinates 
% - ielem_first : global number of the first node in the cell to separate 
% - n           : number of nodes to separate (nodes of the parent cell)
% - xsep        : coordinate of separation 
% - dim         : dimension in which to separate (1 for x, 2 for y)
% - ielem       : current division of elements in the cells

% OUTPUTS 
% - ielem       : nodes of the tree ordered accordingly to the new separation 
% - nsep        : nodes in the lower part of the separation  

nsep = 0;                             % number of nodes that are separated (in the lower part)  
index = ielem_first;   % index of first element of parent cell

% % Check that there's at least one node to separate 
% if n <= 0
%    disp('Error: no nodes to separate') 
% end

% % % For every node of the cell  
% % for node = 1:n 
% %     % If the coord is less than the separation value move the elem at the
% %     % start of the ielem vector 
% %     if x(dim,ielem(index + node -1)) <= xsep      
% %         % Swap the elements 
% %         aux                    = ielem(index);
% %         ielem(index)           = ielem(index + node -1);
% %         ielem(index + node -1) = aux;
% %         % Increase the count 
% %         nsep = nsep + 1;
% %     end
% % end

% MY VERSION 
if n > 0 
    for node = 1:n 
        % If the coord is less than the separation value move the elem at the
        % start of the ielem vector 
        if x(dim,ielem(index + node -1)) <= xsep      
            % Swap the elements 
            if node == 1
                % no need to change it is already at the start
            elseif node == n 
                % move last elem of section to the start of section
                ielem(index:index+n-1) = [ielem(index+n-1);ielem(index:index+n-2)];
            else
                % move the element at the start of section 
                ielem(index:index+n-1) = [ielem(index+node-1);ielem(index:index+node-2);ielem(index+node:index+n-1)]; 
            end
            nsep = nsep + 1;
        end
    end
end

end