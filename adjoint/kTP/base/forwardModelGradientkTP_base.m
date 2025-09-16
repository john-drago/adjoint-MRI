function [ partialfpartialp, varInfo, initSt ] = ...
    forwardModelGradientkTP_base( partialfpartialp, initSt, wv, opt, tt )
% This function will calculate the gradp_fArray at a given time interval.

%% Calculate partialfpartialB function first
% initSt = handlepartialfpartialBRotationMatrix( opt, initSt, wv, Mtt, Mttp1, tt );

%% Initialize varInfo struct
varInfo = initializeVarInfoForwardModelGradient( opt );

%% Get partialB partial parameter for the given time point 

%% RF constant period (also depends on dt)
if any( tt == opt.RF_idx, "all" ) % RF period

    RF_idx = find( tt == opt.RF_idx );
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

elseif any( tt == opt.RF_Slew_i_idx, "all" ) % RF slew up period

    RF_Slew_i_idx = find( tt == opt.RF_Slew_i_idx );
    ktp_b_idx = ( RF_Slew_i_idx - 1 ) * opt.numXYCoils + (1:opt.numXYCoils).';

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

elseif any( tt == opt.RF_Slew_f_idx, "all" ) % RF slew down period
    

    RF_Slew_f_idx = find( tt == opt.RF_Slew_f_idx );
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
elseif any( tt == opt.blip_idx, "all" ) % blip period

    blip_idx = find( tt == opt.blip_idx );

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
    
end

end
