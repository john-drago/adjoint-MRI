function [ h, gradxT_h, gradp_h ] = runMagnitudeLeastSquaresCost( Mfinal, ~, opt )
% This function will run the magnitude least squares optimization for the
% universal or tailored pulse optimization

errNRMSEVec = zeros( opt.numSubj, 1, 'like', Mfinal );

if nargout > 1
    gradxT_h = zeros( opt.cum_numPosSubj(end), 3, 'like', Mfinal );
end

for ss = 1:opt.numSubj

    if ss == 1
        subjIdxs = 1 : opt.cum_numPosSubj( ss );
    else
        subjIdxs = ( opt.cum_numPosSubj( ss-1 ) + 1 ) : opt.cum_numPosSubj( ss );
    end

    Mxy = sqrt( ( Mfinal( subjIdxs, 1 ) ).^2 + ( Mfinal( subjIdxs, 2 ) ).^2 );
    Mxytarg = sqrt( ( opt.Mtarg( subjIdxs, 1 ) ).^2 + ( opt.Mtarg( subjIdxs, 2 ) ).^2 );
    MxytargMLS = norm( Mxytarg, 2 );
    errxy = Mxy - Mxytarg;
    errMLS = norm( errxy, 2 );

    errNRMSEVec( ss ) = errMLS / MxytargMLS;

    %% Derive gradients
    if nargout > 1
        
        errxydbMxy = errxy ./ Mxy;
        errxydbMxy( isnan(errxydbMxy) ) = 0;
        errxydbMxy( isinf(errxydbMxy) ) = 0;

        errxydbMxydberrMLS = errxydbMxy / errMLS;

        gradxT_h( subjIdxs, 1 ) = ( 1/( MxytargMLS * opt.numSubj ) * errxydbMxydberrMLS ) .* Mfinal( subjIdxs, 1 );
        gradxT_h( subjIdxs, 2 ) = ( 1/( MxytargMLS * opt.numSubj ) * errxydbMxydberrMLS ) .* Mfinal( subjIdxs, 2 );
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