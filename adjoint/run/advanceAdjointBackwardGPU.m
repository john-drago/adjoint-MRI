function lambda_nnm1 = advanceAdjointBackwardGPU( Rnn, lambda_nn, gradx_g_nn )
% This function will integrate lambda backwards in time using the previous
% rotation matrices

if isempty( gradx_g_nn )
    
    numPos = size( Rnn, 1 );
    Rnn = reshape( Rnn, [ numPos, 3, 3 ] );
    lambda_nnm1 = reshape( sum( Rnn .* lambda_nn, 2 ), [numPos, 3] );

else
    error( "Haven't implemented nonzero gradp_g integration yet." )
end
end