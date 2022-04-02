% ORIGINAL 

zi = 1;
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

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Compute the influence of the ntylr coeff of the parent cell
% on the lower degree coeff of the child cell
for l = 0:ntylr
    zi = 1;
    % Each influence only the lower degrees
    for m = l:ntylr
        % Sum to the icell local exp coefficient of order m
        % the influence of the local coeff of order k+m of the parent cell
        b(l+1,icell) = b(l+1,icell) + zi * b(m+1,io);
        zi = zi*z0/(l+1);
    end
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Compute the influence of the ntylr coeff of the parent cell
% on the lower degree coeff of the child cell
for l = 0:ntylr
    % Each influence only the lower degrees
    for m = l:ntylr
        % Sum to the icell local exp coefficient of order m
        % the influence of the local coeff of order k+m of the parent cell
        zi = (z0)^(m-l)/factorial(m-l);
        b(l+1,icell) = b(l+1,icell) + zi * b(m+1,io);
        zi = zi*z0/(l+1);
    end
end