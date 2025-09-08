function [ c, gradc ] = constraintSpokesAvgGlobalSAR( pSc, opt, spokes )
% This function will calculate the average Global SAR for kT-point
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgGlobalSARconstr = opt.avgGlobalSAR_constr;

dutyCycle = opt.dutyCycle;
numSpokesCurr = length( opt.breal_idx ) / opt.numXYCoils;
numXYCoils = opt.numXYCoils;
QGlobal = opt.QGlobal;
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

c_QGlobal = 0;

if nargout > 1
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
end

for nn = 1:numSpokesCurr
    
    spokeIdx = spokes.numSpokes - nn + 1;
    colIdx = numSpokesCurr - nn + 1;

    RF_spoke_base = spokes.RF_spokes{ spokeIdx };
    dtvec_spoke = spokes.dtvec_spokes{ spokeIdx };
    RF_spoke = transpose( bcomp_rshp( :, colIdx ) * RF_spoke_base );

    RF_spoke_pwr = ( ctranspose( RF_spoke ) .* dtvec_spoke ) * RF_spoke;

    c_QGlobal_spoke = real( squeeze( sum( QGlobal .* RF_spoke_pwr, [ 1 2 ] ) ) );

    c_QGlobal = c_QGlobal + c_QGlobal_spoke;

    if nargout > 1

        % Determine contribution to average global SAR from each point
        gradf_zbar_rshp = ( ( (dutyCycle) / pulseLength ) * ( 2 ) )...
            * ( RF_spoke * conj( QGlobal ) );

        gradf_zbar_breal_rshp = ( ( RF_spoke_base .* dtvec_spoke ) *...
            real( gradf_zbar_rshp ) ) .* transpose( breal_sc_rshp( :, colIdx ) );
        gradf_zbar_bimag_rshp = ( ( RF_spoke_base .* dtvec_spoke ) *...
            imag( gradf_zbar_rshp ) ) .* transpose( bimag_sc_rshp( :, colIdx ) );

        % Assign real and imaginary parts of the gradient
        gradc_unsc( breal_idx_rshp( :, colIdx ) ) = gradf_zbar_breal_rshp;
        gradc_unsc( bimag_idx_rshp( :, colIdx ) ) = gradf_zbar_bimag_rshp;

    end

end

c_unsc = ( (dutyCycle) / pulseLength ) * ( c_QGlobal ) - avgGlobalSARconstr;

c = c_unsc / avgGlobalSARconstr;

if nargout > 1
    % assign to gradc
    gradc = gradc_unsc / avgGlobalSARconstr;
    % gradc = gradc_unsc;
end

end