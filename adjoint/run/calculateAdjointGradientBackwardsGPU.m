function [ gradp_g_plus_lambdaTgradp_f, mu_0 ] = ...
    calculateAdjointGradientBackwardsGPU(...
    opt, Marray, Rarray, wv, gradxT_h, g, gradp_g, gradx_g )
%
% This function will perform the backwards time integration of:
%
% dlambda/dt = -( gradx_f lambda + gradx_Jc )
%
% It will also perform the integration of:
%
% \int_0^T ( gradp_Jc + lambda^T gradp_f ) dt
%

% Calculate integral value using right Riemann
gradp_g_plus_lambdaTgradp_f = zeros( opt.numVars, 1, "gpuArray" );

% initialize lambda array
lambda_nn = gpuArray( gradxT_h );
lambda_nnm1 = zeros( size( gradxT_h ), "gpuArray" );

% Struct for initialized variables
partialfpartialBarray = zeros( wv.numPos, opt.numTimePoints, 9, "gpuArray" );
partialfpartialBarray = calcpartialfpartialBRotationArray( partialfpartialBarray, wv, Marray );

% clear Marray;
partialfpartialBarray = permute( partialfpartialBarray, [ 1, 3, 2 ] );
partialfpartialBxarray = partialfpartialBarray( :, 1:3, : );
partialfpartialByarray = partialfpartialBarray( :, 4:6, : );
partialfpartialBzarray = partialfpartialBarray( :, 7:9, : );
clear partialfpartialBarray;
initSt = struct;

% g, gradp_g, gradx_g should already be gpuArrays
% Initialize partialfpartialp array for time point nn
partialfpartialp_zero_nn = zeros( opt.numPos, 3, opt.estMaxActiveVarsTimeStep, "gpuArray" );
lambdaTpartialfpartialp_zeros_nn = zeros( opt.numVars, 1, "gpuArray" );

% set tolerance for dwxy
if ~wv.constantRotatingFrame
    dwxy_tol = 1e-0;
end

% Integrate lambda backwards
for nn = uint32( ( opt.numTimePoints ) : -1 : 1 )

    initSt.partialfpartialBx = partialfpartialBxarray( :, :, nn );
    initSt.partialfpartialBy = partialfpartialByarray( :, :, nn );
    initSt.partialfpartialBz = partialfpartialBzarray( :, :, nn );
    initSt.Mnnp1 = Marray( :, :, (nn + 1));
    initSt.Mnn = Marray( :, :, nn);
    
    partialfpartialp_nn = partialfpartialp_zero_nn;
    [ partialfpartialp_nn, varInfo ] = opt.forwardModelGradientFunction(...
        partialfpartialp_nn, initSt, wv, opt, nn );

    % update df/dp if there was rotating frame change
    if ~wv.constantRotatingFrame && ( nn < opt.numTimePoints )
        dwxydiff = wv.dwxyvec( nn+1 ) - wv.dwxyvec( nn );
        if abs( dwxydiff ) > dwxy_tol
            trot = wv.tvec( nn ) + 0.5 * wv.dtvec( nn );
            partialfpartialp_nn = convertMToRotFrame( partialfpartialp_nn, dwxydiff, trot );
        end
    end

    lambdaTpartialfpartialp_nn = lambdaTpartialfpartialp_zeros_nn;
    lambdaTpartialfpartialp_nn( varInfo.varIdxs( 1:varInfo.varAmtsCtr ) ) = reshape(...
        reshape( lambda_nn, uint32([ 1, opt.numPos*3 ]) ) *...
        reshape( partialfpartialp_nn( :, :, 1:varInfo.varAmtsCtr ), uint32([ opt.numPos*3, varInfo.varAmtsCtr ]) ),...
        uint32([varInfo.varAmtsCtr, 1]) );

    if isempty( g )
        dJdp_contribution_nn = lambdaTpartialfpartialp_nn;
    else
        % Need to determine if next line is accurate given that we are not
        % doing a classic Riemann integration
        dJdp_contribution_nn = gradp_g( :, nn ) + lambdaTpartialfpartialp_nn;
    end
    
    % clear lambdaTpartialfpartialp_nn;

    gradp_g_plus_lambdaTgradp_f = gradp_g_plus_lambdaTgradp_f + ( dJdp_contribution_nn );

    % advance lambda backwards
    % determine if we need to account for rotating frame change in back
    % propagation of lambda_nnp1
    if ~wv.constantRotatingFrame && ( nn < opt.numTimePoints )
        if abs( dwxydiff ) > dwxy_tol
            lambda_nn = convertMFromRotFrame( ...
                lambda_nn, dwxydiff, trot );
        end
    end

    if isempty( g )
        lambda_nnm1 = advanceAdjointBackwardGPU(...
            Rarray( :, :, ( nn ) ), lambda_nn, [] );
    else
        lambda_nnm1 = advanceAdjointBackwardGPU(...
            Rarray( :, :, ( nn ) ), lambda_nn, gradx_g( :, :, (nn) ) );
    end

    lambda_nn = lambda_nnm1;

end

clear Rarray Marray;

if nargout > 1
    mu_0 = -lambda_nnm1;
end

end