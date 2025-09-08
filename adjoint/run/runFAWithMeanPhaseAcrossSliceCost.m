function [ h, gradxT_h, gradp_h ] = runFAWithMeanPhaseAcrossSliceCost( Mfinal, ~, opt )
% This function will run the magnitude least squares optimization for the
% universal or tailored pulse optimization

% errNRMSEVec = zeros( opt.numSubj, 1, 'like', Mfinal );
% 
% if nargout > 1
%     gradxT_Jf = zeros( opt.cum_numPosSubj(end), 3, 'like', Mfinal );
% end

errNRMSEVec = zeros( opt.numSubj, 1 );
if nargout > 1
    gradxT_h = zeros( opt.cum_numPosSubj(end), 3 );
end

zeroTol = eps( 1e1 );

for ss = 1:opt.numSubj

    if ss == 1
        subjIdxs = 1 : opt.cum_numPosSubj( ss );
    else
        subjIdxs = ( opt.cum_numPosSubj( ss-1 ) + 1 ) : opt.cum_numPosSubj( ss );
    end

    numPosSubj = opt.numPosSubj( ss );

    Mx = Mfinal( subjIdxs, 1 );
    My = Mfinal( subjIdxs, 2 );
    Mz = Mfinal( subjIdxs, 3 );
    if opt.useGPU
        Mzbiggerthan1 = Mz > gpuArray( 1 );
        Mzlessthanm1 = Mz < gpuArray( -1 );

        Mz( Mzbiggerthan1 ) = gpuArray( 1 );
        Mz( Mzlessthanm1 ) = gpuArray( -1 );
    else
        Mzbiggerthan1 = Mz > 1;
        Mzlessthanm1 = Mz < -1;

        Mz( Mzbiggerthan1 ) = 1;
        Mz( Mzlessthanm1 ) = -1;
    end
    
    subjSlicesGloIdx = opt.subjSlicesGloIdxs( :, ss );
    subjNSlicesGloIdx = opt.subjNSlicesGloIdxs( :, ss );

    subjSlicesLocIdx = opt.subjSlicesLocIdxs{ ss, 1 };
    subjNSlicesLocIdx = opt.subjNSlicesLocIdxs{ ss, 1 };
    numInSlice = opt.numInSliceSubjs( ss, 1 );
    numOutSlice = opt.numOutSliceSubjs( ss, 1 );

    % uniqueIdxSlice = opt.uniqueIdxs( subjSlicesGloIdx, : );
    % uniqueIdxNSlice = opt.uniqueIdxs( subjNSlicesGloIdx, : );
    uniqueLocIdxSlice = zeros( numInSlice, opt.numUniqueXYPos );
    uniqueLocIdxSliceLinIdx = sub2ind( [ numInSlice, opt.numUniqueXYPos ],...
        transpose( uint32(1:numInSlice) ), opt.uniqueXYPosToXYPos( subjSlicesGloIdx ) );
    uniqueLocIdxSlice( uniqueLocIdxSliceLinIdx ) = ones( numInSlice, 1 );
    uniqueIdxToPos = opt.uniqueXYPosToXYPos( subjSlicesGloIdx );

    Mx_inslice = Mx( subjSlicesLocIdx );
    My_inslice = My( subjSlicesLocIdx );
    Mz_inslice = Mz( subjSlicesLocIdx );
    % Mx_outslice = Mx( subjNSlicesLocIdx );
    % My_outslice = My( subjNSlicesLocIdx );
    Mz_outslice = Mz( subjNSlicesLocIdx );

    pulseFA_inslice = real( acos( Mz_inslice ) );
    pulseFA_outslice = real( acos( Mz_outslice ) );
    theta_outslice = pulseFA_outslice;
    
    Mxyph_inslice = atan2( My_inslice, Mx_inslice );
    expMxyph_inslice = exp( 1j * Mxyph_inslice );

    theta_inslice = pulseFA_inslice .* expMxyph_inslice;

    % Mxymag_outslice = sqrt( Mx_outslice.^2 + My_outslice.^2 );

    % Mxtarg = opt.Mtarg( subjIdxs, 1 );
    % Mytarg = opt.Mtarg( subjIdxs, 2 );
    Mztarg = opt.Mtarg( subjIdxs, 3 );

    targFA = real( acos( Mztarg ) );
    targFA_inslice = targFA( subjSlicesLocIdx );
    targFA_outslice = targFA( subjNSlicesLocIdx );
    thetatarg_outslice = targFA_outslice;

    % Mxtarg_inslice = Mxtarg( subjSlicesLocIdx );
    % Mytarg_inslice = Mytarg( subjSlicesLocIdx );

    % Mxtarg_outslice = Mxtarg( subjNSlicesLocIdx );
    % Mytarg_outslice = Mytarg( subjNSlicesLocIdx );

    % Mxytargmag_inslice = sqrt( Mxtarg_inslice.^2 + Mytarg_inslice.^2 );
    % Mxytargmag_outslice = sqrt( Mxtarg_outslice.^2 + Mytarg_outslice.^2 );

    Mxsum_inslice_xy = transpose( sum( Mx_inslice .* uniqueLocIdxSlice, 1 ) );
    Mysum_inslice_xy = transpose( sum( My_inslice .* uniqueLocIdxSlice, 1 ) );

    numUniqueIdxs = transpose( sum( uniqueLocIdxSlice, 1 ) );
    invNumUniqueIdxs = 1 ./ numUniqueIdxs;

    Mxavg_inslice_xy = invNumUniqueIdxs .* Mxsum_inslice_xy;
    Myavg_inslice_xy = invNumUniqueIdxs .* Mysum_inslice_xy;

    Mxytargph_inslice = atan2( Myavg_inslice_xy, Mxavg_inslice_xy );
    targFAphase_inslice = exp( 1j * Mxytargph_inslice );

    thetatarg_inslice = targFA_inslice .* targFAphase_inslice( uniqueIdxToPos );

    theta = zeros( numPosSubj, 1 );
    thetatarg = zeros( numPosSubj, 1 );
    
    theta( subjNSlicesLocIdx ) = theta_outslice;
    thetatarg( subjNSlicesLocIdx ) = thetatarg_outslice;

    theta( subjSlicesLocIdx ) = theta_inslice;
    thetatarg( subjSlicesLocIdx ) = thetatarg_inslice;

    dtheta = theta - thetatarg;
    dthetaNorm = norm( dtheta, 2 );
    targNorm = norm( targFA, 2 );

    if targNorm == 0
        targNorm = 1;
    end

    errNRMSEVec( ss ) = dthetaNorm / targNorm;

    %% Derive gradients
    if nargout > 1
        
        dcost_derr = conj( dtheta ) / dthetaNorm;

        dcost_derr_inslice = dcost_derr( subjSlicesLocIdx );
        dcost_derr_outslice = dcost_derr( subjNSlicesLocIdx );

        derr_dthetamag_inslice = expMxyph_inslice;
        % derr_dthetamag_outslice = ones( numOutSlice, 1 );

        derr_dthetaph_inslice = 1j * theta_inslice;
        % derr_dthetatargph_inslice = sparse(...
        %     uint32(transpose(1:numInSlice)),...
        %     uniqueIdxToPos,...
        %     -1j * thetatarg_inslice,...
        %     numInSlice, opt.numUniqueXYPos );
        derr_dthetatargph_inslice_t = sparse(...
            uniqueIdxToPos,...
            uint32(transpose(1:numInSlice)),...
            -1j * thetatarg_inslice,...
            opt.numUniqueXYPos, numInSlice );

        dthetamag_dMz_inslice = - ( ones( numInSlice, 1 ) - Mz_inslice.^2 ).^( -0.5 );
        dthetamag_dMz_outslice = - ( ones( numOutSlice, 1 ) - Mz_outslice.^2 ).^( -0.5 );

        one_Mz_inslice = abs( abs( Mz_inslice ) - 1 ) < zeroTol;
        one_Mz_outslice = abs( abs( Mz_outslice ) - 1 ) < zeroTol;

        dthetamag_dMz_inslice( one_Mz_inslice ) = 0;
        dthetamag_dMz_outslice( one_Mz_outslice ) = 0;

        Mxysq_inslice = ( Mx_inslice.^2 + My_inslice.^2 );

        dthetaph_dMx_inslice = ( -My_inslice ./ Mxysq_inslice );
        dthetaph_dMy_inslice = (  Mx_inslice ./ Mxysq_inslice );

        zeroIdx = sqrt( Mxysq_inslice ) < zeroTol;

        dthetaph_dMx_inslice( zeroIdx ) = 0;
        dthetaph_dMy_inslice( zeroIdx ) = 0;

        invNumUniqueIdxsPos = invNumUniqueIdxs( uniqueIdxToPos );
        Mxavg_inslice_xy_pos = Mxavg_inslice_xy( uniqueIdxToPos );
        Myavg_inslice_xy_pos = Myavg_inslice_xy( uniqueIdxToPos );

        Mxyavgsq_inslice = Mxavg_inslice_xy.^2 + Myavg_inslice_xy.^2;
        Mxyavgsq_inslice_pos = Mxyavgsq_inslice( uniqueIdxToPos );
        
        
        % dthetatargph_dMx_inslice = sparse(...
        %     uniqueIdxToPos,...
        %     uint32(transpose(1:numInSlice)),...
        %     ( -invNumUniqueIdxsPos ) .* ( Myavg_inslice_xy_pos ./ Mxyavgsq_inslice_pos ),...
        %     opt.numUniqueXYPos, numInSlice );
        % dthetatargph_dMy_inslice = sparse(...
        %     uniqueIdxToPos,...
        %     uint32(transpose(1:numInSlice)),...
        %     (  invNumUniqueIdxsPos ) .* ( Mxavg_inslice_xy_pos ./ Mxyavgsq_inslice_pos ),...
        %     opt.numUniqueXYPos, numInSlice );

        dthetatargph_dMx_inslice = ( -invNumUniqueIdxsPos ) .* ( Myavg_inslice_xy_pos ./ Mxyavgsq_inslice_pos );
        dthetatargph_dMy_inslice = (  invNumUniqueIdxsPos ) .* ( Mxavg_inslice_xy_pos ./ Mxyavgsq_inslice_pos );

        zeroAvgIdx = sqrt( Mxyavgsq_inslice ) < zeroTol;

        dthetatargph_dMx_inslice( zeroAvgIdx ) = 0;
        dthetatargph_dMy_inslice( zeroAvgIdx ) = 0;

        dcost_dthetamag_inslice = real( dcost_derr_inslice .* derr_dthetamag_inslice );
        % dcost_dthetamag_outslice = real( dcost_derr_outslice .* derr_dthetamag_outslice );
        dcost_dthetamag_outslice = real( dcost_derr_outslice );
        dcost_dthetaph_inslice = real( dcost_derr_inslice .* derr_dthetaph_inslice );

        dcost_dthetatargph_inslice = real( derr_dthetatargph_inslice_t * dcost_derr_inslice );

        dcost_dMx_inslice = dcost_dthetaph_inslice .* dthetaph_dMx_inslice +...
            dcost_dthetatargph_inslice( uniqueIdxToPos ) .* dthetatargph_dMx_inslice;
        dcost_dMy_inslice = dcost_dthetaph_inslice .* dthetaph_dMy_inslice +...
            dcost_dthetatargph_inslice( uniqueIdxToPos ) .* dthetatargph_dMy_inslice;
        dcost_dMz_inslice = dcost_dthetamag_inslice .* dthetamag_dMz_inslice;

        dcost_dMz_outslice = dcost_dthetamag_outslice .* dthetamag_dMz_outslice;

        gradxT_h( subjSlicesGloIdx, 1 ) = ( 1 / ( targNorm * opt.numSubj ) ) * dcost_dMx_inslice;
        gradxT_h( subjSlicesGloIdx, 2 ) = ( 1 / ( targNorm * opt.numSubj ) ) * dcost_dMy_inslice;
        gradxT_h( subjSlicesGloIdx, 3 ) = ( 1 / ( targNorm * opt.numSubj ) ) * dcost_dMz_inslice;

        gradxT_h( subjNSlicesGloIdx, 3 ) = ( 1 / ( targNorm * opt.numSubj ) ) * dcost_dMz_outslice;

    end
end

%% determine gradient of Jf with respect to p
% Jf does not directly depend on any parameters in the p vector
if nargout > 2
    gradp_h = zeros( opt.numVars, 1, 'like', Mfinal );
end

%% Take the mean of the error vector
h = mean( errNRMSEVec );

end