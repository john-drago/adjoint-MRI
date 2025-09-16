function [ sum_g, g, gradx_g, gradp_g ] = runRFPowerChebRunningCost( ~, wv, opt )

% per-time power per quadrature (sum over coils)
brealpwr = sum( wv.breal.^2, 1 );
bimagpwr = sum( wv.bimag.^2, 1 );

powerCoeff = ( wv.dutyCycle * cast( opt.runningCostWeightingCoefficient, 'like', wv.breal ) )...
    / ( cast( 2, 'like', wv.breal ) * wv.Z0 * wv.numTimePoints );

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

    time_idx = repmat( transpose( uint32( 1 + ( 1:opt.numTimePoints ) ) ), [ 1, opt.numCheb_RF, opt.numXYCoils ] );
    brealvar_idx = repmat( reshape( opt.breal_idx, [ 1, opt.numCheb_RF, opt.numXYCoils ] ),...
        [ opt.numTimePoints, 1, 1 ] );
    bimagvar_idx = repmat( reshape( opt.bimag_idx, [ 1, opt.numCheb_RF, opt.numXYCoils ] ),...
        [ opt.numTimePoints, 1, 1 ] );

    gradbreal_g = reshape(...
        transpose( ( cast( 2, 'like', wv.breal ) * powerCoeff * wv.breal ) ),...
        [ opt.numTimePoints, 1, opt.numXYCoils ] )...
        .* wv.Tn( :, uint32( 1:opt.numCheb_RF ) );
    breal_idx = sub2ind( [ opt.numVars, (opt.numTimePoints+1) ],...
        brealvar_idx, time_idx );
    gradp_g( breal_idx ) = gradbreal_g;
    
    gradbimag_g = reshape(...
        transpose( ( cast( 2, 'like', wv.bimag ) * powerCoeff * wv.bimag ) ),...
        [ opt.numTimePoints, 1, opt.numXYCoils ] )...
        .* wv.Tn( :, uint32( 1:opt.numCheb_RF ) );
    bimag_idx = sub2ind( [ opt.numVars, (opt.numTimePoints+1) ],...
        bimagvar_idx, time_idx );
    gradp_g( bimag_idx ) = gradbimag_g;
end

end
