function [ h, gradxT_h, gradp_h ] = runArcsinLeastSquaresCost( Mfinal, ~, opt )
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
    targFA = real( asin(Mxytarg) );
    targFAnorm = norm( targFA, 2 );
    errFA = real( asin(Mxy) ) - targFA;
    errFAnorm = norm( errFA, 2 );

    errNRMSEVec( ss ) = errFAnorm / targFAnorm;

    %% Derive gradients
    if nargout > 1
        
        MxdbMxy = Mfinal( subjIdxs, 1 ) ./ Mxy;
        MydbMxy = Mfinal( subjIdxs, 2 ) ./ Mxy;

        gradxT_h( subjIdxs, 1 ) = ( 1/( targFAnorm * opt.numSubj * errFAnorm ) )...
            .* ( errFA ./ sqrt( ones( opt.numPosSubj( ss ), 1 ) - Mxy.^2 ) .* MxdbMxy );
        gradxT_h( subjIdxs, 2 ) = ( 1/( targFAnorm * opt.numSubj * errFAnorm ) )...
            .* ( errFA ./ sqrt( ones( opt.numPosSubj( ss ), 1 ) - Mxy.^2 ) .* MydbMxy );

        gradxT_h( isnan(gradxT_h) ) = 0;
        gradxT_h( isinf(gradxT_h) ) = 0;
    end
end

%% determine gradient of Jf with respect to p
% h does not directly depend on any parameters in the p vector
if nargout > 2
    gradp_h = zeros( opt.numVars, 1, 'like', Mfinal );
end

%% Take the mean of the error vector
h = mean( errNRMSEVec );

end