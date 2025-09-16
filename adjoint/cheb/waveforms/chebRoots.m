function r = chebRoots( c, dom )
% This function will find the roots of the Chebyshev polynomial within the
% domain
%
% This function draws inspiration from CHEBTECH/ROOTS function in chebfun.
%

[ nc, nw ] = size( c );

if nw > 1
    r = cell( 1, nw );
end

for ww = 1:nw
    % Scale the coefficients
    ci = c( :, ww );
    if nc > 2
        ci = ci / max( abs( chebCoeffToVal( ci ) ) );

        % Calculate Roots
        tmp = getRoots( ci );
    elseif nc == 2
        tmp = -ci(1)/ci(2);
        if (tmp < -1) || (tmp > 1 )
            tmp = [];
        end
    else
        if ci == 0
            tmp = 0;
        else
            tmp = [];
        end
    end

    if nargin > 1
        dd = diff( dom );
        if dd > 0
            tmp = chebInvMap( tmp, dom );
        else
            error( "Needs to be increasing order in domain." )
        end
    end

    if nw > 1
        r{ ww } = tmp;
    else
        r = tmp;
    end
end

% see if we need to map to domain



end
%% Helper Functions
% ----------------------------------------------------------------------- %
function r = getRoots( ci )
tol = 100*eps;

nc = size( ci, 1 );
od = 0.5 * ones( nc-2, 1 );

A = diag( od, 1 ) + diag( od, -1 );
A( 1, 2 ) = 1;
A( end, : ) = A( end, : ) - ( 1/( 2*ci(end) ) ) * transpose( ci( 1:(end-1) ) );

ri = eig( A );

% only save real eigenvalues
r = real( ri( abs( imag( ri ) ) < tol ) );

% only save eigenvalues that are within -1 and 1
r = sort( r( abs( r ) <= ( 1 + tol ) ) );

% determine if there is root at the boundary
if ~isempty(r)
    r( 1 ) = max( r(1), -1 );
    r( end ) = min( r(end), 1 );
else
    r = [];
end

end
% ----------------------------------------------------------------------- %