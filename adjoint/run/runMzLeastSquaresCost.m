function [ h, gradxT_h, gradp_h ] = runMzLeastSquaresCost( Mfinal, ~, opt )
% This function will run the least squares optimization for the
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

    MtargLS = norm( opt.Mtarg( subjIdxs, 3 ), 2 );
    diffz = Mfinal( subjIdxs, 3 ) - opt.Mtarg( subjIdxs, 3 );
    errLS = norm( diffz, 2 );

    errNRMSEVec( ss ) = errLS / MtargLS;

    %% Derive gradients
    if nargout > 1
        gradxT_h( subjIdxs, 3 ) = ( 1/( MtargLS * opt.numSubj * errLS ) ) .* diffz;
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