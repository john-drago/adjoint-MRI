function [ c, gradc ] = constraintSpokesAvgLocalSAR( pSc, opt, spokes )
% This function will calculate the (peak) average Local SAR for optimal control
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgLocalSARconstr = opt.avgLocalSAR_constr;

dutyCycle = opt.dutyCycle;
numSpokesCurr = length( opt.breal_idx ) / opt.numXYCoils;
numXYCoils = opt.numXYCoils;
VOPs = opt.VOPs;
numVOPs = opt.numVOPs;
pulseLength = spokes.pulseLength;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal_idx_rshp = reshape( breal_idx, [ numXYCoils, numSpokesCurr ] );
breal_sc_rshp = reshape( breal_sc, [ numXYCoils, numSpokesCurr ] );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numXYCoils, numSpokesCurr ] );

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag_idx_rshp = reshape( bimag_idx, [ numXYCoils, numSpokesCurr ] );
bimag_sc_rshp = reshape( bimag_sc, [ numXYCoils, numSpokesCurr ] );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numXYCoils, numSpokesCurr ] );

bcomp_rshp = complex( breal_rshp, bimag_rshp );

c_VOP = zeros( numVOPs, 1 );

if nargout > 1
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, numVOPs );
end

for nn = 1:numSpokesCurr
    
    spokeIdx = spokes.numSpokes - nn + 1;
    colIdx = numSpokesCurr - nn + 1;

    RF_spoke_unsc = spokes.RF_spokes{ spokeIdx };
    dtvec_spoke = spokes.dtvec_spokes{ spokeIdx };
    RF_spoke = transpose( bcomp_rshp( :, colIdx ) * RF_spoke_unsc );

    RF_spoke_pwr = ( ctranspose( RF_spoke ) .* dtvec_spoke ) * RF_spoke;

    if ~ismatrix( VOPs )
        c_VOP_spoke = real( squeeze( sum( VOPs .* RF_spoke_pwr, [ 1 2 ] ) ) );
    else
        c_VOP_spoke = real( VOPs .* RF_spoke_pwr );
    end

    c_VOP = c_VOP + c_VOP_spoke;

    if nargout > 1

        % Determine contribution to average local SAR from each point
        if ~ismatrix( VOPs )

            gradf_zbar_rshp = ( ( (dutyCycle) / pulseLength ) * ( 2 ) ) *...
                tensorprod( RF_spoke, conj( VOPs ), 2, 1 );

            gradf_zbar_breal_rshp = reshape(...
                tensorprod( transpose( RF_spoke_unsc .* dtvec_spoke ), real( gradf_zbar_rshp ), 1, 1) .* transpose( breal_sc_rshp( :, colIdx ) ), [ numXYCoils, numVOPs ] );
            gradf_zbar_bimag_rshp = reshape(...
                tensorprod( transpose( RF_spoke_unsc .* dtvec_spoke ), imag( gradf_zbar_rshp ), 1, 1) .* transpose( bimag_sc_rshp( :, colIdx ) ), [ numXYCoils, numVOPs ] );
        else
            gradf_zbar_rshp = ( ( (dutyCycle) / pulseLength ) * ( 2 ) )...
                * ( bcomp_rshp * conj( VOPs ) );

            gradf_zbar_breal_rshp = ( transpose( RF_spoke_unsc .* dtvec_spoke ) *...
                real( gradf_zbar_rshp ) ) .* transpose( breal_sc_rshp( :, colIdx ) );
            gradf_zbar_bimag_rshp = ( transpose( RF_spoke_unsc .* dtvec_spoke ) *...
                imag( gradf_zbar_rshp ) ) .* transpose( bimag_sc_rshp( :, colIdx ) );
        end

        % Assign real and imaginary parts of the gradient
        breal_idx_rshp_rep = reshape( repmat( transpose( breal_idx_rshp( :, colIdx ) ), [ 1, 1, numVOPs ] ), [ numXYCoils, numVOPs ] );
        bimag_idx_rshp_rep = reshape( repmat( transpose( bimag_idx_rshp( :, colIdx ) ), [ 1, 1, numVOPs ] ), [ numXYCoils, numVOPs ] );
        vop_idx_rshp_rep = reshape( repmat( reshape( 1:numVOPs, [ 1, 1, numVOPs  ] ), [ 1, numXYCoils ] ), [ numXYCoils, numVOPs ] );
        lin_breal_idx = sub2ind( size( gradc_unsc ), breal_idx_rshp_rep(:), vop_idx_rshp_rep(:) );
        lin_bimag_idx = sub2ind( size( gradc_unsc ), bimag_idx_rshp_rep(:), vop_idx_rshp_rep(:) );
        gradc_unsc( lin_breal_idx ) = gradf_zbar_breal_rshp( : );
        gradc_unsc( lin_bimag_idx ) = gradf_zbar_bimag_rshp( : );

    end

end

c_unsc = ( (dutyCycle) / pulseLength ) * ( c_VOP ) - avgLocalSARconstr;

c = c_unsc / avgLocalSARconstr;

if nargout > 1
    % assign to gradc
    gradc = gradc_unsc / avgLocalSARconstr;
    % gradc = gradc_unsc;
end

end