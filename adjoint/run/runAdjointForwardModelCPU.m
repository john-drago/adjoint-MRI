function [ Marray, Rarray, wv ] = runAdjointForwardModelCPU( wv, opt, saveIntermediate )
% This function will run the CPU implementation of the adjoint forward
% model. It will produce a cell array with opt.numTimePoints+1 entries.

if nargin < 3
    saveIntermediate = true;
end

% initialize Marray struct
if saveIntermediate
    Marray = repmat( {zeros( size( opt.M0 ) )}, ( opt.numTimePoints + 1 ), 1 );
    Marray{ 1 } = opt.M0;
else
    Mnn = opt.M0;
    Mnnp1 = zeros( size( opt.M0 ) );
end

% Initialize Binit so don't need to reinitialize
Binit = zeros( size( opt.M0 ) );

% initialize sensitivity arrays needed
wv.bzsens = opt.bzsens;
wv.b1preal = opt.b1preal;
wv.b1pimag = opt.b1pimag;
wv.db0 = opt.db0;
wv.pos = opt.pos;

% start in correct rotating frame if necessary
if wv.constantRotatingFrame
    if wv.dwxyvec( 1 ) ~= 0
        wv.db0 = wv.db0 + ( 1 / wv.gyro ) * wv.dwxyvec( 1 );
    end
end

% If going to save the rotation matrices for adjoint computation
if nargout > 1

    % initialize rotation matrices
    Rarray = repmat( {zeros( size( opt.M0, 1 ), 9 )}, opt.numTimePoints, 1 );

    for nn = uint32( 1:opt.numTimePoints )
        Rarray{ nn } = generateRotationArray_nn( Binit, Rarray{ nn }, wv, nn );
        Marray{ (nn+1) } = advanceMForwardCPU( Marray{ nn }, Rarray{ nn }, Marray{ (nn+1) }, nn, opt, wv );
    end

else % Don't need every rotation matrix 

    % initialize rotation matrix
    Rarray = zeros( size( opt.M0, 1 ), 9 );

    if saveIntermediate
        for nn = uint32( 1:opt.numTimePoints )
            Rarray = generateRotationArray_nn( Binit, Rarray, wv, nn );
            [ Marray{ (nn+1) } ] = advanceMForwardCPU( Marray{ nn }, Rarray, Marray{ (nn+1) }, nn, opt, wv );
        end

    else
        for nn = uint32( 1:opt.numTimePoints )
            Rarray = generateRotationArray_nn( Binit, Rarray, wv, nn );
                
            [ Mnnp1 ] = advanceMForwardCPU( Mnn, Rarray, Mnnp1, nn, opt, wv );
            Mnn = Mnnp1;
        end

        Marray = Mnnp1;
        clear Mnnp1 Mnn;
    end
end

if isfield( opt, 'convertMBackToLarmor' ) && opt.convertMBackToLarmor
    if saveIntermediate
        Marray{ end } = convertMFromRotFrame( Marray{ end }, wv.dwxyvec( end ), wv.tvec( end ) + 0.5 * wv.dtvec( end ) );
        for nn = uint32( ( opt.numTimePoints ):-1:1 )
            Marray{ nn } = convertMFromRotFrame( Marray{ nn }, wv.dwxyvec( nn ), wv.tvec( nn ) - 0.5 * wv.dtvec( nn ) );
        end
    else
        Marray = convertMFromRotFrame( Marray{ end }, wv.dwxyvec( end ), wv.tvec( end ) + 0.5 * wv.dtvec( end ) );
    end
end

end

%% Helper Function
% ----------------------------------------------------------------------- %
function [ Mnnp1 ] = advanceMForwardCPU( Mnn, Rnn, Mnnp1, nn, opt, wv )

dwxy_tol = 1e-1;

Mnnp1 = applyRotationCPU( Rnn, Mnn, Mnnp1 );

if ~wv.constantRotatingFrame
    if nn < opt.numTimePoints
        dwxydiff = opt.dwxyvec( nn+1 ) - opt.dwxyvec( nn );
        if abs( dwxydiff ) > dwxy_tol
            Mnnp1 = convertMToRotFrame(...
                Mnnp1, dwxydiff, wv.tvec( nn ) + 0.5 * wv.dtvec( nn ) );
        end
    end
end

end
% ----------------------------------------------------------------------- %