function c = chebCoeffFromFun( fn, dom, tol, nmax )
% This function will take a function handle and a domain and determine the
% coefficients from the chebyshev projection.
%
% Inspiration from CHEBTECH2/POPULATE and CHEBTECH2/REFINE in CHEBFUN.

if nargin < 2
    dom = [ -1 1 ];
end
if nargin < 3
    tol = eps;
end
if nargin < 4
    nmax = 5e3; % maximum number of coefficients
end

n = 17; % initially sample 17 points

while true
    
    xi = chebPts2( n );
    xi = chebInvMap( xi, dom );

    cproj = chebValToCoeff( fn( xi ) );
    
    ncut = determineChebCutoff( cproj, tol );

    if ncut < n
        c = cproj( 1:ncut );
        break;
    else
        if n > nmax
            c = cproj( 1:nmax );
            warning( "Could not resolve function completely with %i cheb coefficients", nmax );
            break;
        end
        n = 2*n - 1;
    end

end

end