function [ c, gradc ] = constraintMPpTxPulseShimTotal( pSc, opt )
% This function will calculate the total shim constraint value and (possibly
% the gradient) based on the design vector, dSc.

shimTotalConstr = opt.shimTotal_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

% Constraint during the blip period
ablip_idx = opt.shim_Blip_idx;
ablip_sc = opt.scVec( ablip_idx );
ablip = pSc( opt.shim_Blip_idx ) .* ablip_sc;
cblip = norm( ablip, 1 ) - shimTotalConstr;

% Constraint during the multiphoton subpulse
are_idx = opt.shimreal_MPSP_idx;
are_sc = opt.scVec( are_idx );
aim_idx = opt.shimimag_MPSP_idx;
aim_sc = opt.scVec( aim_idx );
are = pSc( are_idx ) .* are_sc;
aim = pSc( aim_idx ) .* aim_sc;
amag = sqrt( are.^2 + aim.^2 );

cmpsp = -shimTotalConstr + sum( amag );

c = [ cblip; cmpsp ];
c = c / shimTotalConstr;

if nargout > 1
    
    % Initialize gradient vector
    % gradc = sparse( opt.numVars, 2 );
    gradc = zeros( opt.numVars, 2 );

    % gradient of blip constraints
    gradc( ablip_idx, 1 ) = sign( ablip ) .* ablip_sc;
    
    % gradient of mag constraints
    gradc( are_idx, 2 ) = ( are ./ amag ) .* are_sc;
    gradc( aim_idx, 2 ) = ( aim ./ amag ) .* aim_sc;

    gradc( isnan( gradc ) ) = 0;

    gradc = gradc / shimTotalConstr;

end
end