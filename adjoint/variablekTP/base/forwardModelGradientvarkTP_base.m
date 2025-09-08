function [ partialfpartialp, varInfo, initSt ] = ...
    forwardModelGradientvarkTP_base( partialfpartialp, initSt, wv, opt, nn )
% This function will calculate the gradp_fArray at a given time interval.

%% Calculate partialfpartialB function first
% initSt = handlepartialfpartialBRotationMatrix( opt, initSt, wv, Mnn, Mnnp1, nn );

%% Initialize varInfo struct
varInfo = initializeVarInfoForwardModelGradient( opt );

%% Get timing information
dt_RF_idx = opt.dt_idx( 1:2:end );
dt_blip_idx = opt.dt_idx( 2:2:end );

%% RF constant period (also depends on dt)
if any( nn == opt.RF_idx, "all" ) % RF period

    RF_idx = find( nn == opt.RF_idx );
    ktp_b_idx = ( RF_idx - 1 ) * opt.numXYCoils + (1:opt.numXYCoils).';

    % breal contribution
    breal_idxs = opt.breal_idx( ktp_b_idx );
    varInfo = addVarIdxs( breal_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        initSt.partialfpartialBx .* reshape( wv.b1preal, [ opt.numPos, 1, opt.numXYCoils ] ) +...
        initSt.partialfpartialBy .* reshape( wv.b1pimag, [ opt.numPos, 1, opt.numXYCoils ] ),...
        varInfo );

    % bimag contribution
     bimag_idxs = opt.bimag_idx( ktp_b_idx );
     varInfo = addVarIdxs( bimag_idxs, varInfo );

     partialfpartialp = placepartialfpartialpCalculation(...
         partialfpartialp,...
         initSt.partialfpartialBx .* reshape( ( - wv.b1pimag ), [ opt.numPos, 1, opt.numXYCoils ] ) +...
         initSt.partialfpartialBy .* reshape( wv.b1preal, [ opt.numPos, 1, opt.numXYCoils ] ),...
         varInfo );

    %% Deal with dt
    [ partialfpartialphi, magBnn ] = determine_partialfpartialphi( opt.useGPU, nn, wv, initSt.Mnnp1 );

    varInfo = addVarIdxs( dt_RF_idx( RF_idx ), varInfo );
    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        ( -wv.gyro * magBnn ) .* partialfpartialphi,...
        varInfo );

elseif any( nn == opt.RF_Slew_i_idx, "all" ) % RF slew up period

    RF_Slew_i_idx = find( nn == opt.RF_Slew_i_idx );

    ktp_b_idx = ( RF_Slew_i_idx - 1 ) * wv.numXYCoils + (1:wv.numXYCoils).';
        
    % breal contribution
    breal_idxs = opt.breal_idx( ktp_b_idx );
    varInfo = addVarIdxs( breal_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        0.5 * (initSt.partialfpartialBx .* reshape( wv.b1preal, [ opt.numPos, 1, opt.numXYCoils ] ) +...
        initSt.partialfpartialBy .* reshape( wv.b1pimag, [ opt.numPos, 1, opt.numXYCoils ] ) ),...
        varInfo );

    % bimag contribution
     bimag_idxs = opt.bimag_idx( ktp_b_idx );
     varInfo = addVarIdxs( bimag_idxs, varInfo );

     partialfpartialp = placepartialfpartialpCalculation(...
         partialfpartialp,...
         0.5 * (initSt.partialfpartialBx .* reshape( ( - wv.b1pimag ), [ opt.numPos, 1, opt.numXYCoils ] ) +...
         initSt.partialfpartialBy .* reshape( wv.b1preal, [ opt.numPos, 1, opt.numXYCoils ] ) ),...
         varInfo );

elseif any( nn == opt.RF_Slew_f_idx, "all" ) % RF slew down period

    RF_Slew_f_idx = find( nn == opt.RF_Slew_f_idx );

    ktp_b_idx = ( RF_Slew_f_idx - 1 ) * wv.numXYCoils + (1:wv.numXYCoils).';
        
    % breal contribution
    breal_idxs = opt.breal_idx( ktp_b_idx );
    varInfo = addVarIdxs( breal_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        0.5 * (initSt.partialfpartialBx .* reshape( wv.b1preal, [ opt.numPos, 1, opt.numXYCoils ] ) +...
        initSt.partialfpartialBy .* reshape( wv.b1pimag, [ opt.numPos, 1, opt.numXYCoils ] ) ),...
        varInfo );

    % bimag contribution
     bimag_idxs = opt.bimag_idx( ktp_b_idx );
     varInfo = addVarIdxs( bimag_idxs, varInfo );

     partialfpartialp = placepartialfpartialpCalculation(...
         partialfpartialp,...
         0.5 * (initSt.partialfpartialBx .* reshape( ( - wv.b1pimag ), [ opt.numPos, 1, opt.numXYCoils ] ) +...
         initSt.partialfpartialBy .* reshape( wv.b1preal, [ opt.numPos, 1, opt.numXYCoils ] ) ),...
         varInfo );

%% Deal with Blip period
elseif any( nn == opt.blip_idx, "all" ) % blip period

    blip_idx = find( nn == opt.blip_idx );
    
    % Grad
    ktp_grad_idx = ( blip_idx - 1 ) * 3 + (1:3).';
    grad_idxs = opt.grad_idx( ktp_grad_idx );
    varInfo = addVarIdxs( grad_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        calcpartialfpartialG( 0.5 * wv.pos, initSt.partialfpartialBz, opt ),...
        varInfo );
    
    if opt.numZCoils > 0
        % Shim
        ktp_shim_idx = ( blip_idx - 1 ) * opt.numZCoils + (1:opt.numZCoils).';
        shim_idxs = opt.shim_idx( ktp_shim_idx );
        varInfo = addVarIdxs( shim_idxs, varInfo );

        partialfpartialp = placepartialfpartialpCalculation(...
            partialfpartialp,...
            calcpartialfpartialshim( 0.5 * wv.bzsens, initSt.partialfpartialBz, opt ),...
            varInfo );
    end

    %% Deal with dt

    [ partialfpartialphi, magBnn ] = determine_partialfpartialphi( opt.useGPU, nn, wv, initSt.Mnnp1 );

    varInfo = addVarIdxs( dt_blip_idx( blip_idx ), varInfo );
    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        ( -wv.gyro * magBnn ) .* partialfpartialphi,...
        varInfo );
    

end

end

%% Helper Function
% ----------------------------------------------------------------------- %
function [ partialfpartialphi, magBnn ] = determine_partialfpartialphi( useGPU, nn, wv, Mnnp1 )

Bx = wv.b1preal * wv.breal( :, nn ) - wv.b1pimag * wv.bimag( :, nn );
By = wv.b1pimag * wv.breal( :, nn ) + wv.b1preal * wv.bimag( :, nn );
Bz = wv.bzsens * wv.Shim( :, nn ) + wv.pos * wv.Grad( :, nn ) + wv.db0;

magBnn = sqrt( Bx.^2 + By.^2 + Bz.^2 );

ux = Bx ./ magBnn ;
clear Bx;
uy = By ./ magBnn ;
clear By;
uz = Bz ./ magBnn ;
clear Bz;

if useGPU
    partialfpartialphi = zeros( wv.numPos, 3, "gpuArray" );
end

partialfpartialphi( :, 1 ) = - uz .* Mnnp1( :, 2 ) + uy .* Mnnp1( :, 3 );
partialfpartialphi( :, 2 ) = + uz .* Mnnp1( :, 1 ) - ux .* Mnnp1( :, 3 );
partialfpartialphi( :, 3 ) = - uy .* Mnnp1( :, 1 ) + ux .* Mnnp1( :, 2 );


end
% ----------------------------------------------------------------------- %