function [ c, gradc ] = constraintvarkTPPulseAvgGlobalSAR( pSc, opt )
% This function will calculate the average Global SAR for kT-point
% pulse and (possibly the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

avgGlobalSARconstr = opt.avgGlobalSAR_constr;

dutyCycle = opt.dutyCycle;
num_kTP = opt.num_kTP;
numXYCoils = opt.numXYCoils;
QGlobal = opt.QGlobal;

dt_idx = opt.dt_idx;
dt_sc = opt.scVec( dt_idx );
pulseLength = sum( dt_sc .* pSc( dt_idx ) );

dt_blip_idx = dt_idx( 2:2:end );
dt_blip_sc = opt.scVec( dt_blip_idx );

dt_RF_idx = dt_idx( 1:2:end );
dt_RF_sc = opt.scVec( dt_RF_idx );
dt_RF = ( dt_RF_sc .* pSc( dt_RF_idx ) ).';
dt_RF_noslew = dt_RF - 2 * opt.RFSlewTime;

breal_idx = opt.breal_idx;
breal_sc = opt.scVec( breal_idx );
% breal_idx_rshp = reshape( breal_idx, [ numXYCoils, num_kTP ] );
% breal_sc_rshp = reshape( breal_sc, [ numXYCoils, num_kTP ] );
breal = breal_sc .* pSc( breal_idx );
breal_rshp = reshape( breal, [ numXYCoils, num_kTP ] );
breal_half_rshp = breal_rshp / 2;

bimag_idx = opt.bimag_idx;
bimag_sc = opt.scVec( bimag_idx );
% bimag_idx_rshp = reshape( bimag_idx, [ numXYCoils, num_kTP ] );
% bimag_sc_rshp = reshape( bimag_sc, [ numXYCoils, num_kTP ] );
bimag = bimag_sc .* pSc( bimag_idx );
bimag_rshp = reshape( bimag, [ numXYCoils, num_kTP ] );
bimag_half_rshp = bimag_rshp / 2;

bcomp_rshp = complex( breal_rshp, bimag_rshp ).';
bcomp_rshp_t = transpose( bcomp_rshp );
bcomp_half_rshp = complex( breal_half_rshp, bimag_half_rshp ).';

b_pwr_high = ( bcomp_rshp' .* dt_RF_noslew ) * bcomp_rshp;
b_pwr_half = ( bcomp_half_rshp' .* (2 * opt.RFSlewTime ) ) * bcomp_half_rshp;

b_pwr = b_pwr_high + b_pwr_half;

avgGlobalSAR = real( squeeze( sum( QGlobal .* b_pwr, [ 1 2 ] ) ) );

c_unsc = ( (dutyCycle) ) * ( avgGlobalSAR ) - avgGlobalSARconstr * pulseLength;

% c = c_unsc / avgLocalSARconstr;
c = c_unsc;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );
    
    % Take gradient of f = z^H Q z according to CR calculus:
    % gradf_zbar = Q z, such that: df = 2 * Re( gradf_zbar^T dz )
    % gradf_z = Q z, such that: df = 2 * Re( gradf_z^H dz )
    
    % Determine contribution to average local SAR from each point
    gradf_zbar_b_rshp = ( ( (dutyCycle) ) * ( 2 ) )...
        * ( ( bcomp_rshp .* ( dt_RF_noslew.' + (0.5 * opt.RFSlewTime ) ) ) * conj( QGlobal ) );
    gradf_zbar_b = reshape( transpose( gradf_zbar_b_rshp ), [ num_kTP * numXYCoils ,1 ] );

    % Assign real and imaginary parts of the gradient
    gradc_unsc( breal_idx ) = real( gradf_zbar_b ) .* breal_sc;
    gradc_unsc( bimag_idx ) = imag( gradf_zbar_b ) .* bimag_sc;

    % Determine contribution to average local SAR from dt
    gradc_dt_unsc = dutyCycle * ( diag( real( bcomp_rshp_t' * ( QGlobal * bcomp_rshp_t ) ) ) ...
        .* dt_RF_sc ) -...
        avgGlobalSARconstr * dt_RF_sc;
    gradc_unsc( dt_RF_idx ) = gradc_dt_unsc;

    gradc_unsc( dt_blip_idx, 1) = - avgGlobalSARconstr * dt_blip_sc;
    
    % assign to gradc
    % gradc = gradc_unsc / avgGlobalSARconstr;
    gradc = gradc_unsc;



    
    % % Determine contribution to average local SAR from each kT-point
    % for kk = 1:num_kTP
    %     gradf_zbar = ( ( (dutyCycle) / numTimePoints ) * ( 2 * ( RFIntPtNum + 2 * sum( slewScale.^2 ) ) ) )...
    %         * ( QGlobal * ( bcomp_rshp( kk, : ) ).' );
    % 
    %     % Determine variables indices for pTx channels corresponding to peak
    %     % kT-point
    %     var_kTP_idx = (kk-1) * numXYCoils + (1:numXYCoils).';
    % 
    %     % Assign real and imaginary parts of the gradient
    %     gradc_unsc( breal_ktp_idx( var_kTP_idx ) ) = real( gradf_zbar ) .* breal_ktp_sc( var_kTP_idx );
    %     gradc_unsc( bimag_ktp_idx( var_kTP_idx ) ) = imag( gradf_zbar ) .* bimag_ktp_sc( var_kTP_idx );
    % 
    % end
    % 
    % gradc = gradc_unsc / avgGlobalSARconstr;
end
end