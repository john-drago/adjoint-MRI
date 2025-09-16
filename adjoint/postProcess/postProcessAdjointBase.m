function [ val, opt ] = postProcessAdjointBase( val, opt, pulse, fields )

%% Get optimization variables
val.pOpt = opt.pOpt;
val.pScOpt = opt.pScOpt;

%% Simulate resulting magnetization
% First obtain waveform from the opt struct
wvPlotOpt = val.generatePlotWaveforms( opt.pOpt, opt );

% interpolate waveforms to get validation waveforms at possibly a different
% dt.
wvVal = wvPlotOpt;
wvVal.dt = val.dt;
if isfield( val, 'get_tvec' )
    val = val.get_tvec( val.pOpt, val );
    wvVal.tvec = val.tvec;
    wvVal.dtvec = val.dtvec;
else
    wvVal.tvec = val.tvec;
    wvVal.dtvec = val.dtvec;
end
wvVal.numPos = val.numPos;
wvVal.numVars = val.numVars;
wvVal.numTimePoints = length( wvVal.tvec );
wvVal.dwxyvec = val.dwxyvec;

wvVal.brealphasor = interpolateWaveforms( wvPlotOpt.tvec, real(wvPlotOpt.RFphasor), wvVal.tvec );
wvVal.bimagphasor = interpolateWaveforms( wvPlotOpt.tvec, imag(wvPlotOpt.RFphasor), wvVal.tvec );
wvVal.breal = interpolateWaveforms( wvPlotOpt.tvec, real( wvPlotOpt.RF ), wvVal.tvec );
wvVal.bimag = interpolateWaveforms( wvPlotOpt.tvec, imag( wvPlotOpt.RF ), wvVal.tvec );
wvVal.Grad = interpolateWaveforms( wvPlotOpt.tvec, wvPlotOpt.Grad, wvVal.tvec );
wvVal.Shim = interpolateWaveforms( wvPlotOpt.tvec, wvPlotOpt.Shim, wvVal.tvec );

wvVal = rmfield( wvVal, "RF" );
wvVal = rmfield( wvVal, "RFphasor" );
wvVal = rmfield( wvVal, "Freq" );

% Run forward model integration
simTimeTic = tic;
saveIntermediate = false;
Mopt = runAdjointForwardModelCPU( wvVal, val, saveIntermediate );
val.simTime = toc( simTimeTic );
% Mopt = Marray{ end };

clear wvPlotOpt;

%% Calculate Power and Energy Outputs
[ val ] = processAdjointPower( val, wvVal );

%% Calculate Slew and Acceleration
[ val ] = processAdjointSlewAccel( val, wvVal );

clear wvVal;

%% Calculate Flip Angle Distribution
FAopt = real( acosd( Mopt( :, 3 ) ) );

%% Calculate Birdcage Hardpulse Performance
if ( val.numXYCoils == 1 ) && val.performBCHPcomp
    tBCHP = val.pulseLength;
    val.tBCHP = tBCHP;

    [ val.BCHP_BCmag, FABodyBCHP, MBodyBCHP, val.BCHP_RFPower ]...
        = BCHPforEachMemberDB( val );
end

%% Regnerate Mxy and Mz in the domain after Bloch Sim
arraySz = [ length( val.x ), length( val.y ), length( val.z ) ];

val.opt_roi_array = false( [ arraySz, val.numSubj ] );
val.M_targ_array = zeros( [ arraySz, 3, val.numSubj ] );
val.FA_targ_array = zeros( [ arraySz, val.numSubj ] );

val.M_opt_array = zeros( [ arraySz, 3, val.numSubj ] );
val.FA_opt_array = zeros( [ arraySz, val.numSubj ] );

if ( val.numXYCoils == 1 ) && val.performBCHPcomp
    val.M_BCHP_array = zeros( [ arraySz, 3, val.numSubj ] );
    val.FA_BCHP_array = zeros( [ arraySz, val.numSubj ] );
end

for nn = 1:val.numSubj

    if nn == 1
        NonOptIdxs = 1:( val.cum_numPosSubj(nn));
    else
        NonOptIdxs = ( val.cum_numPosSubj(nn-1)+1 ):( val.cum_numPosSubj(nn) );
    end
    
    % repopulate body positions
    [ ri, ci, di ] = ind2sub( arraySz, val.idx(NonOptIdxs) );

    li = sub2ind( [ arraySz, val.numSubj ], ri, ci, di, nn*ones(size(ri)) );
    xi = sub2ind( [ arraySz, 3, val.numSubj ], ri, ci, di, 1*ones(size(ri)), nn*ones(size(ri)) );
    yi = sub2ind( [ arraySz, 3, val.numSubj ], ri, ci, di, 2*ones(size(ri)), nn*ones(size(ri)) );
    zi = sub2ind( [ arraySz, 3, val.numSubj ], ri, ci, di, 3*ones(size(ri)), nn*ones(size(ri)) );

    val.M_targ_array( xi ) = val.Mtarg( NonOptIdxs, 1 );
    val.M_targ_array( yi ) = val.Mtarg( NonOptIdxs, 2 );
    val.M_targ_array( zi ) = val.Mtarg( NonOptIdxs, 3 );
    val.FA_targ_array( li ) = real( acosd( val.Mtarg( NonOptIdxs, 3 ) ) );

    val.M_opt_array( xi ) = Mopt( NonOptIdxs, 1 );
    val.M_opt_array( yi ) = Mopt( NonOptIdxs, 2 );
    val.M_opt_array( zi ) = Mopt( NonOptIdxs, 3 );
    val.FA_opt_array( li ) = FAopt( NonOptIdxs );

    if ( val.numXYCoils == 1 ) && val.performBCHPcomp
        val.M_BCHP_array( xi ) = MBodyBCHP( NonOptIdxs, 1 );
        val.M_BCHP_array( yi ) = MBodyBCHP( NonOptIdxs, 2 );
        val.M_BCHP_array( zi ) = MBodyBCHP( NonOptIdxs, 3 );
        val.FA_BCHP_array( li ) = FABodyBCHP( NonOptIdxs );
    end

    % repopulate brain positions
    optIdxs = NonOptIdxs( val.optroi(NonOptIdxs) );
    [ ri, ci, di ] = ind2sub( arraySz, val.idx( optIdxs  ) );
    li = sub2ind( [ arraySz, val.numSubj ], ri, ci, di, nn*ones(size(ri)) );

    val.opt_roi_array( li ) = true( size( li ));

end

%% Generate vecs for metrics
M_targ_vec_opt_db = cell( val.numSubj, 1 );
FA_targ_vec_opt_db = cell( val.numSubj, 1 );

M_vec_opt_db = cell( val.numSubj, 1 );
FA_vec_opt_db = cell( val.numSubj, 1 );

if ( val.numXYCoils == 1 ) && val.performBCHPcomp
    M_vec_BCHP_db = cell( val.numSubj, 1 );
    FA_vec_BCHP_db = cell( val.numSubj, 1 );
end

numPosSubj = zeros( val.numSubj, 1 );

for nn = 1:val.numSubj

    % roi
    roiSubj = val.opt_roi_array( :, :, :, nn );
    numPts = sum( roiSubj, "all" );
    numPosSubj( nn ) = numPts;
    falseROIsubj = false( size( roiSubj ) );
    xROI = cat( 4, roiSubj, falseROIsubj, falseROIsubj );
    yROI = cat( 4, falseROIsubj, roiSubj, falseROIsubj );
    zROI = cat( 4, falseROIsubj, falseROIsubj, roiSubj );

    % Mtarg vec
    M_targ_vec_opt_db{ nn } = zeros( numPts, 3 );
    Mtargsubj_opt = val.M_targ_array( :, :, :, :, nn );
    M_targ_vec_opt_db{ nn }( :, 1 ) = Mtargsubj_opt( xROI );
    M_targ_vec_opt_db{ nn }( :, 2 ) = Mtargsubj_opt( yROI );
    M_targ_vec_opt_db{ nn }( :, 3 ) = Mtargsubj_opt( zROI );

    % FAtarg vec
    FA_targsubj_opt = val.FA_targ_array( :, :, :, nn );
    FA_targ_vec_opt_db{ nn } = FA_targsubj_opt( roiSubj );
    
    % M vec
    M_vec_opt_db{ nn } = zeros( numPts, 3 );
    Msubj_opt = val.M_opt_array( :, :, :, :, nn );
    M_vec_opt_db{ nn }( :, 1 ) = Msubj_opt( xROI );
    M_vec_opt_db{ nn }( :, 2 ) = Msubj_opt( yROI );
    M_vec_opt_db{ nn }( :, 3 ) = Msubj_opt( zROI );

    % FA vec
    FAsubj_opt = val.FA_opt_array( :, :, :, nn );
    FA_vec_opt_db{ nn } = FAsubj_opt( roiSubj );
    
    if ( val.numXYCoils == 1 ) && val.performBCHPcomp
        % M BCHP vec
        M_vec_BCHP_db{ nn } = zeros( numPts, 3 );
        Msubj_BCHP = val.M_BCHP_array( :, :, :, :, nn );
        M_vec_BCHP_db{ nn }( :, 1 ) = Msubj_BCHP( xROI );
        M_vec_BCHP_db{ nn }( :, 2 ) = Msubj_BCHP( yROI );
        M_vec_BCHP_db{ nn }( :, 3 ) = Msubj_BCHP( zROI );

        % FA BCHP vec
        FAsubj_BCHP = val.FA_BCHP_array( :, :, :, nn );
        FA_vec_BCHP_db{ nn } = FAsubj_BCHP( roiSubj );
    end

end

%% Concatenate the arrays
Mtarg = cell2mat( M_targ_vec_opt_db(:) );
FAtarg = cell2mat( FA_targ_vec_opt_db(:) );

FAopt = cell2mat( FA_vec_opt_db(:) );
Mopt = cell2mat( M_vec_opt_db(:) );

if ( val.numXYCoils == 1 ) && val.performBCHPcomp
    FABCHP = cell2mat( FA_vec_BCHP_db(:) );
    MBCHP = cell2mat( M_vec_BCHP_db(:) );
end

cum_numPosSubj = cumsum( numPosSubj );

%% Calculate Metrics
% Calculate optimization metrics
if pulse.sliceSelective
    [ val ] = calc2DMetrics( val, pulse, fields );
else
    iden = "";
    [ val ] = calc3DMetrics( val, Mopt, FAopt, Mtarg, FAtarg, iden, cum_numPosSubj );

    if ( val.numXYCoils == 1 ) && val.performBCHPcomp
        % Calculate Birdcage metrics
        iden = "BCHP";
        [ val ] = calc3DMetrics( val, MBCHP, FABCHP, Mtarg, FAtarg, iden, cum_numPosSubj );
    end
end

%% Determine which constraints are active
tolConstraint = 1e-2;
[ opt.AbActive, opt.AbeqActive, opt.lbubActive, opt.nonlconIneqActive, opt.nonlconEqActive ] =...
    determineActiveConstraints(...
    opt.pScOpt, opt, tolConstraint );

val.AbActive = opt.AbActive;
val.AbeqActive = opt.AbeqActive;
val.lbubActive = opt.lbubActive;
val.nonlconIneqActive = opt.nonlconIneqActive;
val.nonlconEqActive = opt.nonlconEqActive;

%% Output Constraints
opt = determinePercentConstraints( val, opt, pulse );

%% Generate plotting arrays
% First obtain waveform from the opt struct
wvVal = val.generatePlotWaveforms( opt.pOpt, opt );

wvVal.numPos = val.numPos;
wvVal.numVars = val.numVars;
wvVal.numTimePoints = length( wvVal.tvec );

val.wv = wvVal;

%% Add values to val structs
% val.Mtarg = Mtarg;
% val.FAtarg = FAtarg;
% 
% val.FAopt = FAopt;
% val.Mopt = Mopt;
% 
% if ( val.numXYCoils == 1 ) && val.performBCHPcomp
%     val.FABCHP = FABCHP;
%     val.MBCHP = MBCHP;
% end

end