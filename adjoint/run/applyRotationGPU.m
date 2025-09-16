function Mnnp1 = applyRotationGPU( Rnn, Mnn )
% this function does fiber matrix multiplication, which is fast on a GPU.
% Compare to applyRotationCPU.

numPos = size( Rnn, 1 );
Rnn = reshape( Rnn, numPos, 3 , 3 );
Mnn = reshape( Mnn, [numPos, 1, 3] );
Mnnp1 = sum( Rnn .* Mnn, 3 );

end