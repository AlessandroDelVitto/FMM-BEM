% ORIGINAL 
% Local coeff of order zero from moment of order zero 

b(1,icell) = b(1,icell) - log(z0) * a(1,jcell);

% First term
zo = 1;

% Compute the influence of every moment of jcell on the
% proper local coeff of icell
for m = 1: nexp+ntylr
    % Divide by the distance between the centers 
    zo = zo/z0;
    % Compute the coeff span that is influenced by Mm
    kmin = max(0,m-nexp);
    kmax = min(m,ntylr);
    % Compute sign of kmin term
    sgn = (-1)^kmin;
    % For every term from kmin to kmax 
    for k = kmin+1:kmax+1
        % M2L formula 
        b(k,icell) = b(k,icell) + sgn*zo*a(m-k,jcell); 
        % Change sign 
        sgn = -sgn;
    end
    % Mltiply for the factorial factor
    zo = zo*m;
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Local coeff of order zero from moment of order zero 
b(1,icell) = b(1,icell) - log(z0) * a(1,jcell);

for l = 0:ntylr
    % Compute sign of kmin term
    % sgn = (-1)^l/(2*pi);
    sgn = (-1)^l;
    for k = 1:nexp
        % M2L formula
        zo = factorial(k+l-1)/z0^(l+k);
        b(l+1,icell) = b(l+1,icell) + sgn*zo*a(k+1,jcell); 
    end
end