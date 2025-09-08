function [ c, gradc ] = constraintMPpTxPulseFixedPhaseShimTotal( pSc, opt )
% This function will calculate the total shim constraint value and (possibly
% the gradient) based on the design vector, dSc.

shimTotalConstr = opt.shimTotal_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

% Constraint during the multiphoton subpulse
amag_idx = opt.shimmag_MPSP_idx;
amag_sc = opt.scVec( amag_idx );

amag = pSc( amag_idx ) .* amag_sc;

cmpsp = -1 + norm( amag, 1 )/shimTotalConstr;

c = cmpsp;

if nargout > 1
    
    % Initialize gradient vector
    gradc = zeros( opt.numVars, 1 );
    
    % gradient of mag constraints
    gradc( amag_idx, 1 ) = ( sign( amag ) .* amag_sc ) / shimTotalConstr;

end
end