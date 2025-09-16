function [ c, gradc ] = constraintvarkTPPulseAvgLocalSAR( pSc, opt )
% This function will calculate the (peak) average Local SAR for kT-point
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgLocalSARconstr = opt.avgLocalSAR_constr;

dutyCycle = opt.dutyCycle;
num_kTP = opt.num_kTP;
numXYCoils = opt.numXYCoils;
VOPs = opt.VOPs;
numVOPs = opt.numVOPs;

dt_idx = opt.dt_idx;
dt_sc = opt.scVec( dt_idx );
pulseLength = sum( dt_sc .* pSc( dt_idx ) );

dt_blip_idx = dt_idx( 2:2:end );
dt_blip_sc = opt.scVec( dt_blip_idx );

dt_RF_idx = dt_idx( 1:2:end );
dt_RF_sc = opt.scVec( dt_RF_idx );
dt_RF = transpose( dt_RF_sc .* pSc( dt_RF_idx ) );
dt_RF_noslew = dt_RF - 2 * opt.RFSlewTime;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
breal_idx_rshp = reshape( breal_idx, [ numXYCoils, num_kTP ] );
breal_sc_rshp = reshape( breal_sc, [ numXYCoils, num_kTP ] );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numXYCoils, num_kTP ] );
% breal_half_rshp = breal_rshp / 2;

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
bimag_idx_rshp = reshape( bimag_idx, [ numXYCoils, num_kTP ] );
bimag_sc_rshp = reshape( bimag_sc, [ numXYCoils, num_kTP ] );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numXYCoils, num_kTP ] );
% bimag_half_rshp = bimag_rshp / 2;

bcomp_rshp_t = complex( breal_rshp, bimag_rshp );
bcomp_rshp = transpose( bcomp_rshp_t );
% bcomp_half_rshp_t = complex( breal_half_rshp, bimag_half_rshp );
% bcomp_half_rshp = transpose( bcomp_half_rshp_t );

% b_pwr_high = ( ctranspose( bcomp_rshp ) .* dt_RF_noslew ) * bcomp_rshp;
% b_pwr_half = ( ctranspose( bcomp_half_rshp ) .* (2 * opt.RFSlewTime ) ) * bcomp_half_rshp;
% b_pwr = b_pwr_high + b_pwr_half;
% avgLocalSAR = real( squeeze( sum( VOPs .* b_pwr, [ 1 2 ] ) ) );

if ~ismatrix( VOPs )
    
    localSARArray_unsc = squeeze( real(...
        sum( conj(bcomp_rshp_t) .* ...
        permute(...
        tensorprod( VOPs, bcomp_rshp_t, 2, 1 ),...
        [ 1, 3, 2 ] ), 1 ) ) );
    if num_kTP > 1
        localSARArray_unsc = transpose( localSARArray_unsc ); 
    end

    localSARArray = localSARArray_unsc .* (dt_RF_noslew + 2 * opt.RFSlewTime * ( 0.5 )^2 );
    avgLocalSAR = sum( localSARArray, 2 );
    
else
    localSARArray_unsc = real( conj(bcomp_rshp_t) .* ( VOPs * bcomp_rshp_t ) );
    localSARArray = localSARArray_unsc .* (dt_RF_noslew + 2 * opt.RFSlewTime * ( 0.5 )^2 );
    avgLocalSAR = sum( localSARArray, 2 );
end

c_unsc = ( (dutyCycle) ) * ( avgLocalSAR ) - avgLocalSARconstr * pulseLength;

% c = c_unsc / avgLocalSARconstr;
c = c_unsc;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, numVOPs );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )

    % Determine contribution to average local SAR from each point
    if ~ismatrix( VOPs )

        gradf_zbar_rshp = ( (dutyCycle) * ( 2 ) ) *...
            tensorprod( ( bcomp_rshp .* ( transpose( dt_RF_noslew ) + 2 * ( opt.RFSlewTime ) * 0.5^2 ) ), conj( VOPs ), 2, 1 );
        gradf_zbar_rshp = permute( gradf_zbar_rshp, [ 2, 1, 3 ] );
    else
        gradf_zbar_rshp = ( ( dutyCycle ) * ( 2 ) )...
            * ( ( bcomp_rshp .* ( transpose( dt_RF_noslew ) + 2 * ( opt.RFSlewTime ) * 0.5^2 ) ) * conj( VOPs ) );
        gradf_zbar_rshp = reshape( gradf_zbar_rshp, [ 1, num_kTP ] );
    end
    
    gradf_zbar_breal_rshp = real( gradf_zbar_rshp ).* breal_sc_rshp;
    gradf_zbar_bimag_rshp = imag( gradf_zbar_rshp ).* bimag_sc_rshp;

    % Assign real and imaginary parts of the gradient
    breal_idx_rshp_rep = reshape( repmat( breal_idx_rshp, [ 1, 1, numVOPs ] ), [ num_kTP*numXYCoils, numVOPs ] );
    bimag_idx_rshp_rep = reshape( repmat( bimag_idx_rshp, [ 1, 1, numVOPs ] ), [ num_kTP*numXYCoils, numVOPs ] );
    vop_b_idx_rshp_rep = reshape( repmat( reshape( 1:numVOPs, [ 1, 1, numVOPs  ] ), [ num_kTP, numXYCoils ] ), [ num_kTP*numXYCoils, numVOPs ] );
    lin_breal_idx = sub2ind( size( gradc_unsc ), breal_idx_rshp_rep(:), vop_b_idx_rshp_rep(:) );
    lin_bimag_idx = sub2ind( size( gradc_unsc ), bimag_idx_rshp_rep(:), vop_b_idx_rshp_rep(:) );
    gradc_unsc( lin_breal_idx ) = gradf_zbar_breal_rshp( : );
    gradc_unsc( lin_bimag_idx ) = gradf_zbar_bimag_rshp( : );

    % Determine contribution to average local SAR from dt
    gradc_dt = localSARArray_unsc .* ( dutyCycle * transpose( dt_RF_sc ) )...
        - avgLocalSARconstr * transpose( dt_RF_sc );
    
    dt_RF_idx_rep = repmat( transpose( dt_RF_idx ), [ numVOPs, 1 ] );
    vop_rfdt_idx_rshp_rep = repmat( reshape( 1:numVOPs, [ numVOPs, 1  ] ), [ 1, num_kTP ] );
    lin_rfdt_idx = sub2ind( size( gradc_unsc ), dt_RF_idx_rep(:), vop_rfdt_idx_rshp_rep(:) );
    gradc_unsc( lin_rfdt_idx ) = gradc_dt( : );

    gradc_unsc( dt_blip_idx, :) = repmat( ( -avgLocalSARconstr * dt_blip_sc ), [ 1, numVOPs ] );
    
    % assign to gradc
    % gradc = gradc_unsc / avgLocalSARconstr;
    gradc = gradc_unsc;

end
end