function [ sum_g, g, gradx_g, gradp_g ] = runRFPowerPWCRunningCost( ~, wv, opt )

% per-time power per quadrature (sum over coils)
brealpwr = sum( wv.breal.^2, 1 );
bimagpwr = sum( wv.bimag.^2, 1 );

powerCoeff = ( wv.dutyCycle * cast( opt.runningCostWeightingCoefficient, 'like', wv.Z0 ) )...
    / ( cast( 2, 'like', wv.Z0 ) * wv.Z0 * wv.numTimePoints );

g = powerCoeff * ( brealpwr + bimagpwr );

g = [ zeros( 1, 'like', g ), g ];

sum_g = sum( g );

if nargout > 2
    % if opt.useGPU
    %     gradx_g = zeros( [ opt.numPos, 3, (opt.numTimePoints+1) ], "gpuArray" );
    % else
    %     gradx_g = repmat( {zeros( size( opt.M0 ) )}, ( opt.numTimePoints + 1 ), 1 );
    % end
    gradx_g = [];

    % deal with gradp_g
    if opt.useGPU
        gradp_g = zeros( [ opt.numVars, (opt.numTimePoints+1) ], "gpuArray" );
    else
        gradp_g = zeros( [ opt.numVars, (opt.numTimePoints+1) ] );
    end

    time_idx = repmat( transpose( uint32( 1 + ( 1:opt.numTimePoints ) ) ), [ 1, opt.numXYCoils ] );

    breal_idx = sub2ind( [ opt.numVars, (opt.numTimePoints+1) ],...
        reshape( opt.breal_idx, [ opt.numTimePoints, opt.numXYCoils ]),...
        time_idx );
    gradbreal_g = ( cast( 2, 'like', wv.breal ) * powerCoeff ) * transpose( wv.breal );
    gradp_g( breal_idx ) = gradbreal_g;

    bimag_idx = sub2ind( [ opt.numVars, (opt.numTimePoints+1) ],...
        reshape( opt.bimag_idx, [ opt.numTimePoints, opt.numXYCoils ]),...
        time_idx );
    gradbimag_g = ( cast( 2, 'like', wv.bimag ) * powerCoeff ) * transpose( wv.bimag );
    gradp_g( bimag_idx ) = gradbimag_g;

end

end
