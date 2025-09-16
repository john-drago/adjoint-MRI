function [ oc, opt, costFn ] = prepareTrackConvergence( oc, opt, costFn )

if isfield( oc, 'trackConvergence' ) && oc.trackConvergence

    if isfield( oc, 'spacingTrackConvergence' )
        opt.spacingTrackConvergence = oc.spacingTrackConvergence;
    else
        defaultSpacingTrackConvergence = 1;
        oc.spacingTrackConvergence = defaultSpacingTrackConvergence;
        opt.spacingTrackConvergence = defaultSpacingTrackConvergence;
    end

    if isfield( oc, 'constrTolSave' )
        opt.constrTolSave = oc.constrTolSave;
    else
        defaultConstrTolSave = 5e-2;
        oc.constrTolSave = defaultConstrTolSave;
        opt.constrTolSave = defaultConstrTolSave;
    end

    const = struct;
    const.spacingTrackConvergence = opt.spacingTrackConvergence;
    const.constrTolSave = opt.constrTolSave;
    const.A = opt.A;
    const.b = opt.b;
    const.Aeq = opt.Aeq;
    const.beq = opt.beq;
    const.nonlcon = opt.nonlcon;

    costFn = @( pSc ) costFnTrackConvergenceWrapper( pSc, costFn, const, const.spacingTrackConvergence, const.constrTolSave );
    oc.trackConvergence = true;
    opt.trackConvergence = true;
else
    costFn = @( pSc ) costFn( pSc );
    oc.trackConvergence = false;
    opt.trackConvergence = false;
end

end
