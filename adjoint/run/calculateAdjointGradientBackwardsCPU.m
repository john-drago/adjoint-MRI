function [ gradp_g_plus_lambdaTgradp_f, mu_0 ] = ...
    calculateAdjointGradientBackwardsCPU(...
    opt, Marray, Rarray, wv, gradxT_h, g, gradp_g, gradx_g )


gradp_g_plus_lambdaTgradp_f = zeros( opt.numVars, 1 );

% initialize lambda array
% lambdaArray = repmat( {zeros( size( opt.M0 ) )}, ( opt.numTimePoints + 1 ), 1 );
% lambdaArray{ ( opt.numTimePoints + 1 ) } = gradxT_Jf;

if isempty( gradx_g )
    lambda_nn = gradxT_h;
else
    lambda_nn = gradxT_h + gradx_g{ end };
end
lambda_nnm1 = zeros( size( gradxT_h ) );

% Struct for initialized variables
initSt = struct;

partialfpartialBarray = zeros( wv.numPos, opt.numTimePoints, 9 );
partialfpartialBarray = calcpartialfpartialBRotationArray(...
    partialfpartialBarray, wv, cell2mat( reshape( Marray, [ 1, 1, wv.numTimePoints+1  ] ) ) );

partialfpartialBarray = permute( partialfpartialBarray, [ 1, 3, 2 ] );
partialfpartialBxarray = partialfpartialBarray( :, 1:3, : );
partialfpartialByarray = partialfpartialBarray( :, 4:6, : );
partialfpartialBzarray = partialfpartialBarray( :, 7:9, : );


% Initialize partialfpartialp array for time point nn
partialfpartialp_zero_nn = zeros( opt.numPos, 3, opt.estMaxActiveVarsTimeStep );
lambdaTpartialfpartialp_zeros_nn = zeros( opt.numVars, 1 );

% set tolerance for dwxy
if ~wv.constantRotatingFrame
    dwxy_tol = 1e-0;
end

% Integrate lambda backwards
for nn = uint32( ( opt.numTimePoints ) : -1 : 1 )

    initSt.partialfpartialBx = partialfpartialBxarray( :, :, nn );
    initSt.partialfpartialBy = partialfpartialByarray( :, :, nn );
    initSt.partialfpartialBz = partialfpartialBzarray( :, :, nn );
    initSt.Mnnp1 = Marray{ nn + 1 };
    initSt.Mnn = Marray{ nn };

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
    
    if ~isempty( g ) && ~isempty( gradp_g )
        dJdp_contribution_nn = gradp_g( :, nn ) + lambdaTpartialfpartialp_nn;
    else
        dJdp_contribution_nn = lambdaTpartialfpartialp_nn;
    end
    
    % Need to take care of gradp_g^(N)
    % Would also need to take care of gradp_h here
    if ~isempty( gradp_g ) && ( nn == uint32( opt.numTimePoints ) )
        dJdp_contribution_nn = dJdp_contribution_nn + gradp_g( :, (nn+1) );
    end

    gradp_g_plus_lambdaTgradp_f = gradp_g_plus_lambdaTgradp_f + ( dJdp_contribution_nn );
    
    % advance lambda backwards
    if ~wv.constantRotatingFrame && ( nn < opt.numTimePoints )
        if abs( dwxydiff ) > dwxy_tol
            lambda_nn = convertMFromRotFrame( ...
                lambda_nn, dwxydiff, trot );
        end
    end

    if ~isempty( g ) && ~isempty( gradx_g )
        lambda_nnm1 = advanceAdjointBackwardCPU(...
            Rarray{ ( nn ) }, lambda_nn, lambda_nnm1, gradx_g{ (nn) } );
    else
        lambda_nnm1 = advanceAdjointBackwardCPU(...
            Rarray{ ( nn ) }, lambda_nn, lambda_nnm1, [] );
    end

    lambda_nn = lambda_nnm1;

end

clear Rarray Marray;

if nargout > 1
    mu_0 = -lambda_nnm1;
end

end