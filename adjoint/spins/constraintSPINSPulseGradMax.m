function [ c, gradc ] = constraintSPINSPulseGradMax( pSc, opt )
% This function will calculate the max gradient value and its derivative for SPINS pulses.

gradMaxConstr = opt.gradMax_constr;

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

% initSlewTime = opt.initSlewTime;
% spins_int_num = opt.spins_int_num;
% spins_int_i = opt.spins_int_i;
% spins_int_f = opt.spins_int_f;
Gxyz_vec = Gxyz_fn( 0, s );

c_unsc = abs( Gxyz_vec ) - gradMaxConstr;

c = c_unsc / gradMaxConstr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 3 );

    param_idx = opt.param_idx;

    % scVec array
    scVec_vec = opt.scVec( param_idx );
    
    dGxyzdparam_vec = transpose( squeeze( dGxyzdparam_fn( 0, s ) ) );
    
    gradc_unsc( param_idx, : ) =  sign( transpose( Gxyz_vec ) ) .* ( dGxyzdparam_vec .* scVec_vec );

    gradc = gradc_unsc / gradMaxConstr;

end

end