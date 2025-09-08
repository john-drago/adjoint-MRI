function [ opt, oc ] = pSc0Process( opt, oc, ensureFeasibleStart )
if isfield( oc, "pSc0" ) || isfield( opt, "pSc0" )
    if isfield( oc, "pSc0" ) && ~isfield( opt, "pSc0" )
        opt.pSc0 = oc.pSc0;
    elseif isfield( opt, "pSc0" ) && isfield( oc, "pSc0" )
        oc.pSc0 = opt.pSc0;
    elseif isfield( opt, "pSc0" ) && ~isfield( oc, "pSc0" )
        oc.pSc0 = opt.pSc0;
    end

    opt.p0 = opt.scVec .*  opt.pSc0;

    if ensureFeasibleStart
        st = struct;
        if isfield( oc, 'ensureFeasibleStartRelLineSrchBnd' )
            st.relLineSrchBnd = oc.ensureFeasibleStartRelLineSrchBnd;
        end
        if isfield( oc, 'ensureFeasibleStartScaleAtEnd' )
            st.scaleAtEnd = oc.ensureFeasibleStartScaleAtEnd;
        end
        
        [ opt.pSc0, ~ ] = ensureFeasibleStartFMIN( opt.pSc0, opt, st );
        
        opt.p0 = opt.scVec .*  opt.pSc0;
    end

else
    % warning("Initializing pSc0 with zero vector of size:\t[%i,1]",...
    %     opt.numVars);
    % opt.pSc0 = zeros( opt.numVars, 1 );

    pScMax = 1e-4;
    warning("Initializing pSc0 with randn (mean: %g) vector of size:\t[%i,1]",...
        pScMax, opt.numVars);

    opt.pSc0 = pScMax * randn( opt.numVars, 1 );
    if ensureFeasibleStart
        [ opt.pSc0, ~ ] = ensureFeasibleStartFMIN( opt.pSc0, opt );
        opt.p0 = opt.scVec .*  opt.pSc0;
    end
end
end