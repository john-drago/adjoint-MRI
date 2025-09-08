function cutoff = determineChebCutoff( c, tol )
% This function will determine the 
%
% Inspiration from CHEBTECH2/STANDARDCHOP in CHEBFUN.

%% Set tolerance
if nargin < 2
    tol = eps;
end

%% Get number of initial coefficients 
n = size( c, 1 );
cutoff = n;

%% Get a monotonically non-increasing vector
envelope = cummax( abs( c ), 'reverse' );
if envelope(1) == 0
    cutoff = 1;
    return
else
    envelope = envelope / envelope(1);
end

%% Determine where the plateau occurs; taken from CHEBTECH2/standardCHOP
for jj = 2:n
    jj2 = round(1.25*jj + 5); 
    if ( jj2 > n )
        % there is no plateau: exit
        return
    end      
    e1 = envelope(jj);
    e2 = envelope(jj2);
    r = 3*(1 - log(e1)/log(tol));
    plateau = (e1 == 0) | (e2/e1 > r);
    if ( plateau )
        % a plateau has been found: go to Step 3
        plateauPoint = jj - 1;
        break
    end
end

%% Adapt cutoff; taken from CHEBTECH2/standardCHOP

if ( envelope(plateauPoint) == 0 )
    cutoff = plateauPoint;
else
    jj3 = sum(envelope >= tol^(7/6));
    if ( jj3 < jj2 )
        jj2 = jj3 + 1;
        envelope(jj2) = tol^(7/6);
    end
    cc = log10(envelope( 1:jj2 ));
    cc = cc(:);
    cc = cc + linspace(0, (-1/3)*log10(tol), jj2 )';
    [~, d] = min(cc);
    cutoff = max(d - 1, 1);
end

end