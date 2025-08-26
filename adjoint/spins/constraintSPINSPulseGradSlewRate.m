function [ c, gradc ] = constraintSPINSPulseGradSlewRate( pSc, opt )
% This function will calculate the max gradient value and its derivative for SPINS pulses.

gradSlewRateConstr = opt.gradSlewRate_constr;

% Assume pSC is scaled to be between -1 and 1
pSc = pSc( : );

% Get SPINS parameters
kmax = opt.min_kmax + opt.scVec( opt.kmax_idx ) * ( pSc( opt.kmax_idx ) + 1 );
a = opt.min_a + opt.scVec( opt.a_idx ) * ( pSc( opt.a_idx ) + 1 );
b = opt.min_b + opt.scVec( opt.b_idx ) * ( pSc( opt.b_idx ) + 1 );
u = opt.min_u + opt.scVec( opt.u_idx ) * ( pSc( opt.u_idx ) + 1 );
v = opt.min_v + opt.scVec( opt.v_idx ) * ( pSc( opt.v_idx ) + 1 );

% Make struct
s = struct;
s.kmax = kmax;
s.a = a;
s.b = b;
s.u = u;
s.v = v;
s.T = opt.Tspins;
s.gyro = opt.gyro;

initSlewTime = opt.initSlewTime;
spins_int_num = opt.spins_int_num;
spins_int_i = opt.spins_int_i;
spins_int_f = opt.spins_int_f;
finalSlewTime = opt.finalSlewTime;

dGxyzdt_rshp = zeros( 3, spins_int_num + 2 );
dGxyzdt_rshp( :, (1 + 1 ):( 1 + spins_int_num ) ) =...
    dGxyzdt_fn( ( opt.tvec( spins_int_i:spins_int_f ) - initSlewTime ), s );
dGxyzdt_rshp( :, 1 ) = Gxyz_fn( 0, s ) / initSlewTime;
dGxyzdt_rshp( :, ( 2 + spins_int_num ) ) = Gxyz_fn( s.T, s ) / finalSlewTime;
dGxyzdt_vec = reshape( transpose( dGxyzdt_rshp ), [ ( 3*( spins_int_num + 2 ) ), 1 ] );

c_val = abs( dGxyzdt_vec );

c_unsc = c_val - gradSlewRateConstr;
c = c_unsc / gradSlewRateConstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, ( 3 * spins_int_num + 2*3 ) );

    param_idx = opt.param_idx;

    % scVec array
    scVec_arr = zeros( 3, ( spins_int_num + 2 ), length( param_idx ) );
    scVec_arr( :, :, 1 ) = repmat( opt.scVec( param_idx( 1 ) ), [ 3, ( spins_int_num + 2 ) ] );
    scVec_arr( :, :, 2 ) = repmat( opt.scVec( param_idx( 2 ) ), [ 3, ( spins_int_num + 2 ) ] );
    scVec_arr( :, :, 3 ) = repmat( opt.scVec( param_idx( 3 ) ), [ 3, ( spins_int_num + 2 ) ] );
    scVec_arr( :, :, 4 ) = repmat( opt.scVec( param_idx( 4 ) ), [ 3, ( spins_int_num + 2 ) ] );
    scVec_arr( :, :, 5 ) = repmat( opt.scVec( param_idx( 5 ) ), [ 3, ( spins_int_num + 2 ) ] );
    
    ddGxyzdtdparam_arr = zeros( 3, ( spins_int_num + 2 ), length( param_idx )  );
    ddGxyzdtdparam_arr( :, (1 + 1 ):( 1 + spins_int_num ), : ) =...
        ddGxyzdtdparam_fn( ( opt.tvec( spins_int_i:spins_int_f ) - initSlewTime ), s );
    ddGxyzdtdparam_arr( :, 1, : ) = dGxyzdparam_fn( 0, s ) / initSlewTime;
    ddGxyzdtdparam_arr( :, 2 + spins_int_num, : ) = dGxyzdparam_fn( s.T, s ) / finalSlewTime;

    sign_arr = sign( dGxyzdt_rshp );

    grad_unsc_rshp = sign_arr .* ddGxyzdtdparam_arr .* scVec_arr;
    grad_unsc_vec = transpose( reshape( permute( grad_unsc_rshp, [ 2 1 3 ] ), [ (3*(spins_int_num + 2)), length( param_idx ) ] ) );
    
    gradc_unsc( param_idx, : ) = grad_unsc_vec;

    gradc = gradc_unsc / gradSlewRateConstr;

end

end