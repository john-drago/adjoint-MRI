function [ Marray, Rarray, wv ] = runAdjointForwardModelGPU( wv, opt, saveIntermediate )
% This function will run the GPU implementation of the adjoint forward
% model. It will produce a multidimensional array with opt.numTimePoints+1
% entries (last dimension).

if nargin < 3
    saveIntermediate = true;
end

% initialize sensitivity arrays needed
% wv arrays
wv.b1preal = gpuArray( opt.gpu_b1preal );
wv.b1pimag = gpuArray( opt.gpu_b1pimag );
wv.pos = gpuArray( opt.gpu_pos );
wv.bzsens = gpuArray( opt.gpu_bzsens );
wv.db0 = gpuArray( opt.gpu_db0 );

wv.numZCoils = gpuArray( opt.gpu_numZCoils );
wv.numXYCoils = gpuArray( opt.gpu_numXYCoils );
wv.numPos = gpuArray( opt.gpu_numPos );
wv.gyro = gpuArray( opt.gpu_gyro );

wv.tvec = gpuArray( wv.tvec );
wv.dtvec = gpuArray( wv.dtvec );
wv.dwxyvec = gpuArray( wv.dwxyvec );
wv.breal = gpuArray( wv.breal );
wv.bimag = gpuArray( wv.bimag );
wv.Grad = gpuArray( wv.Grad );
wv.Shim = gpuArray( wv.Shim );

% initialize Marray struct
if saveIntermediate
    Marray = zeros( [ size( opt.M0 ), ( opt.numTimePoints + 1 )  ], "gpuArray" );
    Marray( :, :, 1 ) = gpuArray( opt.gpu_M0 );
else
    Mnn = gpuArray( opt.gpu_M0 );
    Mnnp1 = zeros( size( opt.M0 ), "gpuArray" );
end

% start in correct rotating frame if necessary
if wv.constantRotatingFrame
    if opt.dwxyvec( 1 ) ~= 0
        wv.db0 = wv.db0 + ( 1 / wv.gyro ) * wv.dwxyvec( 1 );
    end
end

% Initialize Rarray struct
Rarray = zeros( [ size( opt.M0, 1 ), opt.numTimePoints , 9 ], "gpuArray" );

% Generate rotation matrices for each time point
Rarray = generateRotationArrayGPU( Rarray, wv );
Rarray = permute( Rarray, [ 1 3 2 ] );

% Run forward model
if saveIntermediate
    for nn = uint32( 1:opt.numTimePoints )
        Marray( :, :, (nn+1) ) = advanceMForwardGPU( Marray( :, :, nn ), Rarray( :, :, nn ), nn, opt, wv );
    end
else
    for nn = uint32( 1:opt.numTimePoints )
        Mnnp1 = advanceMForwardGPU( Mnn, Rarray( :, :, nn ), nn, opt, wv );
        Mnn = Mnnp1;
    end
    Marray = Mnnp1;

    clear Mnnp1 Mnn;
end

if isfield( opt, 'convertMBackToLarmor' ) && opt.convertMBackToLarmor
    if saveIntermediate
        Marray( :, :, end ) = convertMFromRotFrame( Marray( :, :, end ), wv.dwxyvec( end ), wv.tvec( end ) + 0.5 * wv.dtvec( end ) );
        for nn = uint32( ( opt.numTimePoints ):-1:1 )
            Marray( :, :, (nn) ) = convertMFromRotFrame( Marray( :, :, (nn) ), wv.dwxyvec( nn ), wv.tvec( nn ) - 0.5 * wv.dtvec( nn ) );
        end
    else
        Marray = convertMFromRotFrame( Marray( :, :, end ), wv.dwxyvec( end ), wv.tvec( end ) + 0.5 * wv.dtvec( end ) );
    end
end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function Mnnp1 = advanceMForwardGPU( Mnn, Rnn, nn, opt, wv )

dwxy_tol = 1e-1;

Mnnp1 = applyRotationGPU( Rnn, Mnn );

if ~wv.constantRotatingFrame
    if nn < opt.numTimePoints
        dwxydiff = wv.dwxyvec( nn+1 ) - wv.dwxyvec( nn );
        if abs( dwxydiff ) > dwxy_tol
            Mnnp1 = convertMToRotFrame(...
                Mnnp1, dwxydiff, wv.tvec( nn ) + 0.5 * wv.dtvec( nn ) );
        end
    end
end

end
% ----------------------------------------------------------------------- %