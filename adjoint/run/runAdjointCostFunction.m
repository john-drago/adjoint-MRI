function [ J, gradpSc_J ] = runAdjointCostFunction( pSc, opt )
% This function will execute the adjoint cost function and calculate the
% gradient if requested

%% Organize parameters
% opt.generateWaveforms typically expects the design vector to be scaled
% between -1 and 1.

scVec = opt.scVec;
p = pSc( : ) .* scVec;

%% Get GPU device information
% if opt.useGPU
%     gpuD = gpuDevice();
% end

%% Generate waveforms
% Generate the waveforms to perform the time domain integration
wv = opt.generateWaveforms( p, opt );

%% Run forward model integration
if nargout > 1

    % Different functions for GPU or CPU implementation
    if opt.useGPU
        [ Marray, Rarray, wv ] = runAdjointForwardModelGPU( wv, opt );
    else
        [ Marray, Rarray, wv ] = runAdjointForwardModelCPU( wv, opt );
    end
else
    if opt.useGPU
        Marray = runAdjointForwardModelGPU( wv, opt );
    else
        Marray = runAdjointForwardModelCPU( wv, opt );
    end
end

%% Calculate terminal cost
if ~isempty( opt.terminalCostFunction )
    % Determine what is the format of the Marray
    if opt.useGPU
        Mfinal = Marray( :, :, (opt.numTimePoints + 1) );
    else
        Mfinal = Marray{ (opt.numTimePoints + 1) };
    end
    
    % Get gradients if necessary for adjoint computation
    if nargout > 1
        [ h, gradxT_h, gradp_h ] = opt.terminalCostFunction( Mfinal, wv, opt );
    else
        h = opt.terminalCostFunction( Mfinal, wv, opt );
    end
else
    h = 0;
    gradxT_h = [];
    gradp_h = [];
end

%% Calculate running cost
if ~isempty( opt.runningCostFunction )
    if nargout > 1
        [ sum_g, g, gradx_g, gradp_g ] = opt.runningCostFunction( Marray, wv, opt );
    else
        sum_g = opt.runningCostFunction( Marray, wv, opt );
    end
else
    sum_g = 0;
    g = [];
    gradp_g = [];
    gradx_g = [];
end

%% Calculate cost function
J = h + sum_g;

if opt.useGPU
    J = gather( J );
end

%% Calculate dJdp (gradient)
% if gradient is requested to be calculated
if nargout > 1

    if opt.useGPU
        [ opt, wv ] = opt.gpuArrayAdjointFunction( opt, wv );
        gradp_g_plus_lambdaTgradp_f = calculateAdjointGradientBackwardsGPU( ...
            opt, Marray, Rarray, wv, gradxT_h, g, gradp_g, gradx_g );
    else
        gradp_g_plus_lambdaTgradp_f = calculateAdjointGradientBackwardsCPU( ...
            opt, Marray, Rarray, wv, gradxT_h, g, gradp_g, gradx_g );
    end

    clear Marray Rarray wv gradxT_h g gradp_g gradx_g;

    gradp_J = gradp_h + gradp_g_plus_lambdaTgradp_f;

    if opt.useGPU
        gradp_J = gather( gradp_J );
    end

    % Scale gradp_J so that it is correct gradient for scaled variables
    gradpSc_J = gradp_J .* scVec;

    gradTol = 1e-10;
    gradpSc_J( abs( gradpSc_J ) < gradTol ) = 0;

    if all( gradpSc_J == 0 )
        dxFD = 1e-6;
        dJdp = calcJacobianFiniteDifference( @(pSc) opt.runCostFunction( pSc(:), opt ), pSc, dxFD );
        gradpSc_J = transpose( dJdp );
    end

end

%% Reset GPU
% if opt.useGPU
%     reset( gpuD );
%     wait( gpuD );
% end

end