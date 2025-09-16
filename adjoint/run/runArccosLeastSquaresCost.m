function [ h, gradxT_h, gradp_h ] = runArccosLeastSquaresCost( Mfinal, ~, opt )
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

    Mztarg = opt.Mtarg( subjIdxs, 3 );
    targFA = real( acos(Mztarg) );
    targFAnorm = norm( targFA, 2 );
    errFA = real(acos(Mz)) - targFA;
    errFAnorm = norm( errFA, 2 );

    errNRMSEVec( ss ) = errFAnorm / targFAnorm;

    %% Derive gradients
    if nargout > 1
        gradxT_h( subjIdxs, 3 ) = ( -1./( targFAnorm * opt.numSubj * errFAnorm ) )...
            .* ( errFA ./ sqrt( ones( opt.numPosSubj( ss ), 1 ) - Mz.^2 ) );

        gradxT_h( isnan(gradxT_h) ) = 0;
        gradxT_h( isinf(gradxT_h) ) = 0;
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