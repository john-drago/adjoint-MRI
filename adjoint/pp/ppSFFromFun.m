function c = ppSFFromFun( wvform, tvec, projSt, tol )

%% Initialize dt tol
dttol = min( diff( tvec ) ) / 10;

wvformsampfunglobal = makeWaveformSampleFunction( tvec, wvform );

%% Initialize arrays that will be used to solve the least-squares projection problem
A = [];
b = [];

%% Iterate over the pulse periods
numPPPeriods = projSt.numPPPeriods;
PPStartStop = projSt.PPStartStop;

PPVarIdxs = projSt.PPVarIdxs;
numVarsPerChannel = projSt.numVarsPerChannel;
numPPShape = projSt.numPPShape;
PPVarSFIdxs = projSt.PPVarSFIdxs;
shapeFnChebCoeffs = projSt.shapeFnChebCoeffs;

wvformbnds = zeros( numPPPeriods, 2 );
for pp = 1:numPPPeriods

    %% Determine bounds of the time domain and the waveform values at the bounds
    tdom_pp = PPStartStop( pp, : );

    if ( pp == 1 ) && ( pp == numPPPeriods )
        wvformbnds( pp, : ) = zeros( 1, 2 );
    elseif pp == 1

        wvformbnds( pp, 1 ) = 0;

        if abs( PPStartStop( pp, 2 ) - PPStartStop( pp+1, 1 ) ) < dttol
            wvformbnds( pp, 2 ) = wvformsampfunglobal( PPStartStop( pp, 2 ) );
        else
            wvformbnds( pp, 2 ) = 0;
        end

    elseif pp < numPPPeriods

        if abs( PPStartStop( pp, 2 ) - PPStartStop( pp+1, 1 ) ) < dttol
            wvformbnds( pp, 2 ) = wvformsampfunglobal( PPStartStop( pp, 2 ) );
        else
            wvformbnds( pp, 2 ) = 0;
        end

    elseif pp == numPPPeriods

        wvformbnds( pp, 2 ) = 0;
    end
    
    %% need to determine tvec indices that correspond to this period
    PPidxs_pp = generatePPtvecIdxs( tdom_pp, tvec, dttol );

    %% create function of waveform over the first domain
    wvform_pp_init = wvform( PPidxs_pp( 1, 1 ):PPidxs_pp( 1, 2 ) ); 
    tvec_pp_init = tvec( PPidxs_pp( 1, 1 ):PPidxs_pp( 1, 2 ) );

    % check bounds
    if abs( tdom_pp( 1 ) - tvec_pp_init( 1 ) ) < dttol
        tvec_pp_init( 1 ) = tdom_pp( 1 );
        wvform_pp_init( 1 ) = wvformbnds( pp, 1 );
    end
    if abs( tdom_pp( 2 ) - tvec_pp_init( end ) ) < dttol
        tvec_pp_init( end ) = tdom_pp( 2 );
        wvform_pp_init( end ) = wvformbnds( pp, 2 );
    end

    tvec_pp = [ tdom_pp( 1 ); tvec_pp_init; tdom_pp( 2 ) ];
    wvform_pp = [ wvformbnds( pp, 1 ); wvform_pp_init; wvformbnds( pp, 2 ) ];

    % check for duplicate points
    [ tvec_pp, tvec_pp_idx ] = unique( tvec_pp, 'sorted' );
    wvform_pp = wvform_pp( tvec_pp_idx );

    wvformsampfun_pp = makeWaveformSampleFunction( tvec_pp, wvform_pp );

    %% Determine chebyshev coefficients necessary to approximate functions in domain
    c_pp = chebCoeffFromFun( wvformsampfun_pp, tdom_pp, tol );
    numc_pp = length( c_pp );

    numc_pp = max( [ numc_pp, double( numPPShape ) ] );

    % determine values at the chebyshev points for Clenshaw-Curtis
    % intergration
    xunsccheb2_pp = chebPts2( numc_pp );
    xcheb2_pp = chebInvMap( xunsccheb2_pp, tdom_pp );
    f_pp = wvformsampfun_pp( xcheb2_pp );

    % determine Clenshaw-Curtis quadrature weights
    wunsc_pp = chebCCQuadWeights( double( numc_pp ) );
    tdomdiff_pp = diff( tdom_pp );
    % scale quadrature weights to account for different domain sizes
    w_pp = ( tdomdiff_pp / 2 ) * wunsc_pp;

    % Calculate contibution to the b vector
    % need to scale the function evaluations by sqrt( diff( tdom_pp ) ) to
    % ensure that we minimize error equally over the different intervals
    % and that we don't artificially weight one interval more than another
    % just because it is longer
    sc_pp = sqrt( w_pp / ( tdomdiff_pp ) );

    b_pp = sc_pp .* f_pp;

    %% Determine contribution to the least squares matrix
    % first determine contribution from variables to shape funs
    varToSF = zeros( numPPShape, numVarsPerChannel );
    varToSFidx = sub2ind( size( varToSF ), PPVarSFIdxs(pp,1):PPVarSFIdxs(pp,2), PPVarIdxs(pp,1):PPVarIdxs(pp,2) );
    varToSF( varToSFidx ) = 1;

    % next determine shapefun value to Tn coefficient value

    shapeFunToTn = shapeFnChebCoeffs;
    TnToIntPts = evalChebClenshaw( xunsccheb2_pp, eye( numPPShape ) );
    varToIntPts = TnToIntPts * shapeFunToTn * varToSF;

    A_pp = sc_pp .* varToIntPts;

    if any( isnan( A_pp ), 'all' )
        error( "NaN in A_pp" );
    end


    %% Add to the arrays for least squares
    A = [ A; A_pp ]; %#ok
    b = [ b; b_pp ]; %#ok

end

%% Solve the linear system
c = linsolve( A, b ); % dense solve

end