function [ partialfpartialp, varInfo, initSt ] = ...
    forwardModelGradientMPpTx_orblipmpconst( partialfpartialp, initSt, wv, opt, tt )
% This function will calculate the gradp_fArray at a given time interval
% for a left Riemann sum.

%% Calculate partialfpartialB function first
% initSt = handlepartialfpartialBRotationMatrix( opt, initSt, wv, Mtt, Mttp1, tt );

%% Initialize varInfo struct
varInfo = initializeVarInfoForwardModelGradient( opt );

%% Get partialB partial parameter for the given time point 
% determine what phase of the pulse we are in

% get time
ti = wv.tvec( tt );
ti_opt = opt.tvec( tt );

%% ORSP
if ( tt >= opt.ORSP_i ) && ( tt <= opt.ORSP_f )
    
   % determine if in slew up period
    if ti_opt < ( opt.tStORSP + opt.RFSlewTime )

        sc = ( (ti - wv.tStORSP ) / wv.RFSlewTime );
        
    % determine if in slew down period
    elseif ti_opt > ( opt.tEndORSP - opt.RFSlewTime )

        sc =( -(ti - wv.tEndORSP ) / wv.RFSlewTime);

    % must be in constant period
    else
        if opt.useGPU
            sc = gpuArray( 1 );
        else
            sc = 1;
        end
    end

    % Calculate contribution from breal_ORSP
    varInfo = addVarIdxs( opt.breal_ORSP_idx, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        sc * (...
        initSt.partialfpartialBx .* reshape( wv.b1preal, [ opt.numPos, 1, opt.numXYCoils ] ) +...
        initSt.partialfpartialBy .* reshape( wv.b1pimag, [ opt.numPos, 1, opt.numXYCoils ] ) ),...
        varInfo );

    % Calculate contribution from bimag_ORSP
    varInfo = addVarIdxs( opt.bimag_ORSP_idx, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        sc * (...
        initSt.partialfpartialBx .* reshape( ( - wv.b1pimag ), [ opt.numPos, 1, opt.numXYCoils ] ) +...
        initSt.partialfpartialBy .* reshape( wv.b1preal, [ opt.numPos, 1, opt.numXYCoils ] ) ),...
        varInfo );
    
end


%% Blip Period
% if in rising Blip Grad Slew
if ( tt >= opt.Blip_i ) && ( tt < opt.Blip_Grad_Slew_i )
    
    grad_idxs = [ opt.Gx_Blip_idx; opt.Gy_Blip_idx; opt.Gz_Blip_idx ];
    varInfo = addVarIdxs( grad_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        calcpartialfpartialG( ( ( ti - wv.tStBlip ) / wv.BlipGradSlewTime ) * wv.pos, initSt.partialfpartialBz, opt ),...
        varInfo );
end

% if in Const Grad
if ( tt >= opt.Blip_Grad_Slew_i ) && ( tt <= opt.Blip_Grad_Slew_f )

    grad_idxs = [ opt.Gx_Blip_idx; opt.Gy_Blip_idx; opt.Gz_Blip_idx ];
    varInfo = addVarIdxs( grad_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        calcpartialfpartialG( wv.pos, initSt.partialfpartialBz, opt ),...
        varInfo );
end

% if in falling Blip Grad Slew
if ( tt <= opt.Blip_f ) && ( tt > opt.Blip_Grad_Slew_f )
    
    grad_idxs = [ opt.Gx_Blip_idx; opt.Gy_Blip_idx; opt.Gz_Blip_idx ];
    varInfo = addVarIdxs( grad_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        calcpartialfpartialG( ( -( ti - wv.tEndBlip ) / wv.BlipGradSlewTime ) * wv.pos, initSt.partialfpartialBz, opt ),...
        varInfo );
end

if opt.numZCoils > 0
    % if in rising Blip Shim Slew
    if ( tt >= opt.Blip_i ) && ( tt < opt.Blip_Shim_Slew_i )
        
        varInfo = addVarIdxs( opt.shim_Blip_idx, varInfo );

        partialfpartialp = placepartialfpartialpCalculation(...
            partialfpartialp,...
            calcpartialfpartialshim( ( (ti - wv.tStBlip )/wv.BlipShimSlewTime ) * wv.bzsens, initSt.partialfpartialBz, opt ),...
            varInfo );

    end

    % if in Const Shim
    if ( tt >= opt.Blip_Shim_Slew_i ) && ( tt <= opt.Blip_Shim_Slew_f )
        varInfo = addVarIdxs( opt.shim_Blip_idx, varInfo );

        partialfpartialp = placepartialfpartialpCalculation(...
            partialfpartialp,...
            calcpartialfpartialshim( wv.bzsens, initSt.partialfpartialBz, opt ),...
            varInfo );
    end

    % if in falling Blip Shim Slew
    if ( tt <= opt.Blip_f ) && ( tt > opt.Blip_Shim_Slew_f )

        varInfo = addVarIdxs( opt.shim_Blip_idx, varInfo );

        partialfpartialp = placepartialfpartialpCalculation(...
            partialfpartialp,...
            calcpartialfpartialshim( ( -(ti - wv.tEndBlip )/wv.BlipShimSlewTime ) * wv.bzsens, initSt.partialfpartialBz, opt ),...
            varInfo );

    end
end


%% MPSP Period
% if in MPSP at all for birdcage coil
if ( tt >= opt.MPSP_i ) && ( tt <= opt.MPSP_f )
    
    % determine if in slew up period
    if ti_opt < ( opt.tStMPSP + opt.RFSlewTime )

        sc = ( (ti - wv.tStMPSP ) / wv.RFSlewTime);
        
    % determine if in slew down period
    elseif ti_opt > ( opt.tEndMPSP - opt.RFSlewTime )
        
        sc = ( -(ti - wv.tEndMPSP ) / wv.RFSlewTime);

    % must be in constant period
    else

        if opt.useGPU
            sc = gpuArray( 1 );
        else
            sc = 1;
        end

    end

    b1prealcos = wv.b1preal * cos( wv.dwxy_mpptx  * ( ti - wv.tStMPSP ) );
    b1pimagsin = wv.b1pimag * sin( wv.dwxy_mpptx  * ( ti - wv.tStMPSP ) );
    b1prealsin = wv.b1preal * sin( wv.dwxy_mpptx  * ( ti - wv.tStMPSP ) );
    b1pimagcos = wv.b1pimag * cos( wv.dwxy_mpptx  * ( ti - wv.tStMPSP ) );

    [ partialfpartialbreal_MPSP_const, partialfpartialbimag_MPSP_const ] = calcpartialfpartialb_MPSP_const(...
        b1prealcos, b1pimagsin, b1prealsin, b1pimagcos, sc, initSt.partialfpartialBx, initSt.partialfpartialBy, opt );

    % Calculate contribution from breal_MPSP
    varInfo = addVarIdxs( opt.breal_MPSP_idx, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        partialfpartialbreal_MPSP_const,...
        varInfo );

    % Calculate contribution from bimag_MPSP
    varInfo = addVarIdxs( opt.bimag_MPSP_idx, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        partialfpartialbimag_MPSP_const,...
        varInfo );

end

% if in rising MPSP Grad Slew
if ( tt >= opt.MPSP_i ) && ( tt < opt.MPSP_Grad_Slew_i )
    
    % Grad
    [ partialfpartialGrise_real, partialfpartialGrise_imag ] = calculatepartialfpartialG_MPSP_rise(...
        initSt.partialfpartialBz, wv, opt, ti );
    
    grad_real_idxs = [ opt.Gxreal_MPSP_idx; opt.Gyreal_MPSP_idx; opt.Gzreal_MPSP_idx ];
    varInfo = addVarIdxs( grad_real_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        partialfpartialGrise_real,...
        varInfo );

    grad_imag_idxs = [ opt.Gximag_MPSP_idx; opt.Gyimag_MPSP_idx; opt.Gzimag_MPSP_idx ];
    varInfo = addVarIdxs( grad_imag_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        partialfpartialGrise_imag,...
        varInfo );

end

% if in Const Grad  
if ( tt >= opt.MPSP_Grad_Slew_i ) && ( tt <= opt.MPSP_Grad_Slew_f )
    
    % Grad
    [ partialfpartialGconst_real, partialfpartialGconst_imag ] = calculatepartialfpartialG_MPSP_const(...
        initSt.partialfpartialBz, wv, opt, ti );
    
    grad_real_idxs = [ opt.Gxreal_MPSP_idx; opt.Gyreal_MPSP_idx; opt.Gzreal_MPSP_idx ];
    varInfo = addVarIdxs( grad_real_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        partialfpartialGconst_real,...
        varInfo );

    grad_imag_idxs = [ opt.Gximag_MPSP_idx; opt.Gyimag_MPSP_idx; opt.Gzimag_MPSP_idx ];
    varInfo = addVarIdxs( grad_imag_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        partialfpartialGconst_imag,...
        varInfo );

end

% if in falling MPSP Grad Slew
if ( tt <= opt.MPSP_f ) && ( tt > opt.MPSP_Grad_Slew_f )
    
    % Grad
    [ partialfpartialGfall_real, partialfpartialGfall_imag ] = calculatepartialfpartialG_MPSP_fall(...
        initSt.partialfpartialBz, wv, opt, ti );
    
    grad_real_idxs = [ opt.Gxreal_MPSP_idx; opt.Gyreal_MPSP_idx; opt.Gzreal_MPSP_idx ];
    varInfo = addVarIdxs( grad_real_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        partialfpartialGfall_real,...
        varInfo );

    grad_imag_idxs = [ opt.Gximag_MPSP_idx; opt.Gyimag_MPSP_idx; opt.Gzimag_MPSP_idx ];
    varInfo = addVarIdxs( grad_imag_idxs, varInfo );

    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        partialfpartialGfall_imag,...
        varInfo );

end

if opt.numZCoils > 0

    % if in rising MPSP Shim Slew
    if ( tt >= opt.MPSP_i ) && ( tt < opt.MPSP_Shim_Slew_i )

        partBpartamag =  ( ( ti - wv.tStMPSP ) / wv.MPSPShimSlewTime ) * cos( wv.wz_mpptx  * wv.MPSPShimSlewTime + wv.shim_ph_MPSP );
        partBpartaph =  ( -( ti - wv.tStMPSP ) / wv.MPSPShimSlewTime ) * (wv.shim_mag_MPSP .* sin( wv.wz_mpptx  * wv.MPSPShimSlewTime + wv.shim_ph_MPSP ) );

        [ partamagpartareal, partaphpartareal, partamagpartaimag, partaphpartaimag ] =...
            calculatepartialshimrealimag_MPSP( wv );
        
        varInfo = addVarIdxs( opt.shimreal_MPSP_idx, varInfo );
        partialfpartialp = placepartialfpartialpCalculation(...
            partialfpartialp,...
            calcpartialfpartialshim(...
            wv.bzsens .* ( ( partBpartamag .* partamagpartareal + partBpartaph .* partaphpartareal ).' ),...
            initSt.partialfpartialBz, opt ),...
            varInfo );

        varInfo = addVarIdxs( opt.shimimag_MPSP_idx, varInfo );
        partialfpartialp = placepartialfpartialpCalculation(...
            partialfpartialp,...
            calcpartialfpartialshim(...
            wv.bzsens .* ( ( partBpartamag .* partamagpartaimag + partBpartaph .* partaphpartaimag ).' ),...
            initSt.partialfpartialBz, opt ),...
            varInfo );
    end

    % if in Const Shim
    if ( tt >= opt.MPSP_Shim_Slew_i ) && ( tt <= opt.MPSP_Shim_Slew_f )

        partBpartamag =   cos( wv.wz_mpptx  * ( ti - wv.tStMPSP ) + wv.shim_ph_MPSP );
        partBpartaph =  ( -wv.shim_mag_MPSP .* sin( wv.wz_mpptx  * ( ti - wv.tStMPSP ) + wv.shim_ph_MPSP ) );

        [ partamagpartareal, partaphpartareal, partamagpartaimag, partaphpartaimag ] =...
            calculatepartialshimrealimag_MPSP( wv );

        varInfo = addVarIdxs( opt.shimreal_MPSP_idx, varInfo );
        partialfpartialp = placepartialfpartialpCalculation(...
            partialfpartialp,...
            calcpartialfpartialshim(...
            wv.bzsens .* ( ( partBpartamag .* partamagpartareal + partBpartaph .* partaphpartareal ).' ),...
            initSt.partialfpartialBz, opt ),...
            varInfo );

        varInfo = addVarIdxs( opt.shimimag_MPSP_idx, varInfo );
        partialfpartialp = placepartialfpartialpCalculation(...
            partialfpartialp,...
            calcpartialfpartialshim(...
            wv.bzsens .* ( ( partBpartamag .* partamagpartaimag + partBpartaph .* partaphpartaimag ).' ),...
            initSt.partialfpartialBz, opt ),...
            varInfo );
    end

    % if in falling MPSP Shim Slew
    if ( tt <= opt.MPSP_f ) && ( tt > opt.MPSP_Shim_Slew_f )

        partBpartamag =  ( -( ti - wv.tEndMPSP ) / wv.MPSPShimSlewTime ) *...
            cos( wv.wz_mpptx  * ( wv.tEndMPSP - opt.tStMPSP - wv.MPSPShimSlewTime ) + wv.shim_ph_MPSP );
        partBpartaph =  ( ( ti - wv.tEndMPSP ) / wv.MPSPShimSlewTime ) *...
            (wv.shim_mag_MPSP .* sin( wv.wz_mpptx  * ( wv.tEndMPSP - wv.tStMPSP - wv.MPSPShimSlewTime ) + wv.shim_ph_MPSP ) );

        [ partamagpartareal, partaphpartareal, partamagpartaimag, partaphpartaimag ] =...
            calculatepartialshimrealimag_MPSP( wv );

        varInfo = addVarIdxs( opt.shimreal_MPSP_idx, varInfo );
        partialfpartialp = placepartialfpartialpCalculation(...
            partialfpartialp,...
            calcpartialfpartialshim(...
            wv.bzsens .* ( ( partBpartamag .* partamagpartareal + partBpartaph .* partaphpartareal ).' ),...
            initSt.partialfpartialBz, opt ),...
            varInfo );

        varInfo = addVarIdxs( opt.shimimag_MPSP_idx, varInfo );
        partialfpartialp = placepartialfpartialpCalculation(...
            partialfpartialp,...
            calcpartialfpartialshim(...
            wv.bzsens .* ( ( partBpartamag .* partamagpartaimag + partBpartaph .* partaphpartaimag ).' ),...
            initSt.partialfpartialBz, opt ),...
            varInfo );
    end
end


end