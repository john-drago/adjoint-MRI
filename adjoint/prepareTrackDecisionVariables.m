function [ oc, opt, costFn ] = prepareTrackDecisionVariables( oc, opt, costFn )

if isfield( oc, 'trackDecisionVariables' ) && oc.trackDecisionVariables

    if isfield( oc, 'spacingTrackDecisionVariables' )
        opt.spacingTrackDecisionVariables = oc.spacingTrackDecisionVariables;
    else
        defaultSpacingTrackDecisionVariables = 10;
        oc.spacingTrackDecisionVariables = defaultSpacingTrackDecisionVariables;
        opt.spacingTrackDecisionVariables = defaultSpacingTrackDecisionVariables;
    end

    if isfield( oc, 'constrTolSave' )
        opt.constrTolSave = oc.constrTolSave;
    else
        defaultConstrTolSave = 5e-2;
        oc.constrTolSave = defaultConstrTolSave;
        opt.constrTolSave = defaultConstrTolSave;
    end

    const.spacingTrackDecisionVariables = opt.spacingTrackDecisionVariables;
    const.constrTolSave = opt.constrTolSave;
    const.A = opt.A;
    const.b = opt.b;
    const.Aeq = opt.Aeq;
    const.beq = opt.beq;
    const.nonlcon = opt.nonlcon;

    costFn = @( pSc ) costFnTrackDecisionVariablesWrapper( pSc, costFn, const, const.spacingTrackDecisionVariables, const.constrTolSave );
    oc.trackDecisionVariables = true;
    opt.trackDecisionVariables = true;
else
    costFn = @( pSc ) costFn( pSc );
    oc.trackDecisionVariables = false;
    opt.trackDecisionVariables = false;
end

end
