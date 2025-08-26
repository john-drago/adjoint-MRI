function note = collateAdjointData( valtrain, valtest, opt, pulse, oc )
% This function will add information to a note that will contain data bout
% the optimization process. This will be used to save a txt file.

%% Initialize Note
note = "ADJOINT PULSE DATA:";
note = [ note; "" ];
note = [ note; opt.timeStamp ];
note = [ note; sprintf( "Iden:\t%s", opt.timeStampIden ) ];
note = [ note; "" ];
note = [ note; "" ];

if isfield( oc, 'note' )
    note = [ note; sprintf( "Note:" ) ];
    note = [ note; oc.note( : ) ];
    note = [ note; "" ];
    note = [ note; "" ];
end

%% Add Path Information
note = [ note; "Path Information:" ];
note = [ note; sprintf( "subjectMaps:\t%s", pulse.subjectMaps ) ];
note = [ note; sprintf( "zCoilPath:\t%s", pulse.zCoilPath ) ];
note = [ note; sprintf( "zCoilShiftPath:\t%s", pulse.zCoilShiftPath ) ];
note = [ note; sprintf( "VOPpath:\t%s", pulse.VOPpath ) ];
note = [ note; "" ];
note = [ note; "" ];

%% Add Subject Information
note = [ note; "Subject Information:" ];
note = [ note; sprintf( "Num Train:\t%i", length( valtrain.subjIden ) ) ];
note = [ note; ...
    strcat( sprintf( "Train Subj:\t[" ), sprintf( " %s ", valtrain.subjIden ), sprintf("]") ) ];
if ~isempty( valtest )
    note = [ note; sprintf( "Num Test:\t%i", length( valtest.subjIden ) ) ];
    note = [ note; ...
        strcat( sprintf( "Test Subj:\t[" ), sprintf( " %s ", valtest.subjIden ), sprintf("]") ) ];
else
    note = [ note; sprintf( "Test Subj:\t%s", "" ) ];
end
note = [ note; "" ];
note = [ note; "" ];

%% Add Pulse Information
note = [ note; "Pulse Information:" ];
note = [ note; sprintf( "pulse.name:\t%s", pulse.name ) ];
note = [ note; sprintf( "pulse.type:\t%s", pulse.type ) ];
note = [ note; sprintf( "pulse.optPopulation:\t%s", pulse.optPopulation ) ];
note = [ note; sprintf( "pulse.excitationType:\t%s", pulse.excitationType ) ];
note = [ note; sprintf( "pulse.terminalCostFunction:\t%s", pulse.terminalCostFunction ) ];
note = [ note; sprintf( "pulse.runningCostFunction:\t%s", pulse.runningCostFunction ) ];

%% Add Pulse Specific Information
note = [ note; sprintf( "pulse.dutyCycle:\t%g", pulse.dutyCycle ) ];
note = [ note; sprintf( "pulse.Z0:\t%g", pulse.Z0 ) ];
note = [ note; sprintf( "pulse.length:\t%g", pulse.length ) ];
if isfield( pulse, 'targFAVal' )
    note = [ note; sprintf( "pulse.targFAVal:\t%g", pulse.targFAVal ) ];
end
note = [ note; sprintf( "opt_dt:\t%g", oc.opt_dt ) ];
note = [ note; sprintf( "val_dt:\t%g", oc.val_dt ) ];
if isfield( opt, 'dfxy' )
    note = [ note; sprintf( "dfxy:\t%g", opt.dfxy ) ];
end
if isfield( opt, 'mdfxy' )
    note = [ note; sprintf( "mdfxy:\t%g", opt.mdfxy ) ];
end
if isfield( opt, 'convertMBackToLarmor' )
    note = [ note; sprintf( "convertMBackToLarmor:\t%g", opt.convertMBackToLarmor ) ];
end
if isfield( opt, 'convertMtargAwayLarmor' )
    note = [ note; sprintf( "convertMtargAwayLarmor:\t%g", opt.convertMtargAwayLarmor ) ];
end
if isfield( opt, 'convertMtargAwayLarmor' )
    note = [ note; sprintf( "convertMtargAwayLarmor:\t%g", opt.convertMtargAwayLarmor ) ];
end
if isfield( opt, 'constantRotatingFrame' )
    note = [ note; sprintf( "constantRotatingFrame:\t%g", opt.constantRotatingFrame ) ];
end

if pulse.sliceSelective
    note = [ note; ...
        strcat( sprintf( "Slice Direction:\t[" ), sprintf( " %g ", pulse.sliceDirection ), sprintf("]") ) ];
    note = [ note; ...
        strcat( sprintf( "Slice Location:\t[" ), sprintf( " %g ", pulse.sliceLocation ), sprintf("]") ) ];
    note = [ note; sprintf( "Slice Thickness:\t%g", pulse.sliceThickness ) ];
    note = [ note; sprintf( "Num Slices:\t%g", pulse.numSlices ) ];
    if isfield( pulse, 'numSpokes' )
        note = [ note; sprintf( "Num Spokes:\t%g", pulse.numSpokes ) ];
    end
    if isfield( pulse, 'centralSpokeTBW' )
        note = [ note; sprintf( "Central Spoke TBW:\t%g", pulse.centralSpokeTBW ) ];
    end
    if isfield( pulse, 'centralSpokeLength' )
        note = [ note; sprintf( "Central Spoke Length:\t%g", pulse.centralSpokeLength ) ];
    end
    if isfield( pulse, 'nonCentralSpokeTBW' )
        note = [ note; sprintf( "Non-Central Spoke TBW:\t%g", pulse.nonCentralSpokeTBW ) ];
    end
    if isfield( pulse, 'nonCentralSpokeLength' )
        note = [ note; sprintf( "Non-Central Spoke Length:\t%g", pulse.nonCentralSpokeLength ) ];
    end
end

if ~isstruct( oc.opt_di ) && ( numel( oc.opt_di ) == 3 )
    note = [ note; ...
        strcat( sprintf( "opt_di:\t[" ), sprintf( " %g ", oc.opt_di ), sprintf("]") ) ];
elseif isstruct( oc.opt_di )
    if isfield( oc.opt_di, 'withinSlice_di' )
        note = [ note; ...
            sprintf( "opt_di.withinSlice_di:\t%g", oc.opt_di.withinSlice_di ) ];
    end
    if isfield( oc.opt_di, 'outsideSlice_di' )
        note = [ note; ...
            sprintf( "opt_di.outsideSlice_di:\t%g", oc.opt_di.outsideSlice_di ) ];
    end
    if isfield( oc.opt_di, 'outsideSlice_di_start' )
        note = [ note; ...
            sprintf( "opt_di.outsideSlice_di_start:\t%g", oc.opt_di.outsideSlice_di_start ) ];
    end
    if isfield( oc.opt_di, 'outsideSlice_di_end' )
        note = [ note; ...
            sprintf( "opt_di.outsideSlice_di_end:\t%g", oc.opt_di.outsideSlice_di_end ) ];
    end
    if isfield( oc.opt_di, 'intraSlice_di' )
        note = [ note; ...
            sprintf( "opt_di.intraSlice_di:\t%g", oc.opt_di.intraSlice_di ) ];
    end
    if isfield( oc.opt_di, 'outsideSliceExt' )
        note = [ note; ...
            sprintf( "opt_di.outsideSliceExt:\t%g", oc.opt_di.outsideSliceExt ) ];
    end
end
if ~isstruct( oc.val_di ) && ( numel( oc.val_di ) == 3 )
    note = [ note; ...
        strcat( sprintf( "val_di:\t[" ), sprintf( " %g ", oc.val_di ), sprintf("]") ) ];
elseif isstruct( oc.val_di )
    if isfield( oc.val_di, 'withinSlice_di' )
        note = [ note; ...
            sprintf( "val_di.withinSlice_di:\t%g", oc.val_di.withinSlice_di ) ];
    end
    if isfield( oc.val_di, 'outsideSlice_di' )
        note = [ note; ...
            sprintf( "val_di.outsideSlice_di:\t%g", oc.val_di.outsideSlice_di ) ];
    end
    if isfield( oc.val_di, 'intraSlice_di' )
        note = [ note; ...
            sprintf( "val_di.intraSlice_di:\t%g", oc.val_di.intraSlice_di ) ];
    end
    if isfield( oc.val_di, 'outsideSliceExt' )
        note = [ note; ...
            sprintf( "val_di.outsideSliceExt:\t%g", oc.val_di.outsideSliceExt ) ];
    end
end

pulseInfoKeys = keys( pulse.pulseinfo );
pulseInfoVals = values( pulse.pulseinfo );
for kk = 1:length( pulseInfoKeys )
    pulseInfoVal = pulse.pulseinfo( pulseInfoKeys( kk ) );
    if ~iscell( pulseInfoVals )
        note = [ note; ...
            sprintf( "%s:\t%s", pulseInfoKeys( kk ), mat2str( pulseInfoVal ) ) ]; %#ok
    else
        note = [ note; ...
            sprintf( "%s:\t%s", pulseInfoKeys( kk ), mat2str( pulseInfoVal{ 1 } ) ) ]; %#ok
    end
end
note = [ note; "" ];
note = [ note; "" ];

%% Add Constraint Information
note = [ note; "Constraint Information:" ];
constraintList = keys( pulse.constraints );
for kk = 1:length( constraintList )
    note = [ note; ...
        sprintf( "%s:\t%g", constraintList( kk ), pulse.constraints( constraintList( kk ) ) ) ]; %#ok
end
note = [ note; "" ];
note = [ note; "" ];

%% Add Optimization Information
note = [ note; "Optimization Information:" ];
note = [ note; sprintf( "oc.saveDir:\t%s", oc.saveDir ) ];
if isfield( oc, 'ensureFeasibleStart' )
    note = [ note; sprintf( "oc.ensureFeasibleStart:\t%g", oc.ensureFeasibleStart ) ];
end
note = [ note; sprintf( "oc.useGPU:\t%g", oc.useGPU ) ];
note = [ note; sprintf( "oc.optType:\t%s", oc.optType ) ];
note = [ note; sprintf( "oc.numWorkers:\t%g", oc.numWorkers ) ];
if matches( oc.optType, "fmin", 'ignorecase', true )
    note = [ note; "" ];
    note = [ note; "fmincon options:" ];
    note = [ note; sprintf( "fminopt.UseParallel:\t%s", string( oc.fminopt.UseParallel ) ) ];
    note = [ note; sprintf( "fminopt.MaxIterations:\t%s", string( oc.fminopt.MaxIterations ) ) ];
    note = [ note; sprintf( "fminopt.MaxFunctionEvaluations:\t%s", string( oc.fminopt.MaxFunctionEvaluations ) ) ];
    note = [ note; sprintf( "fminopt.StepTolerance:\t%s", string( oc.fminopt.StepTolerance ) ) ];
    note = [ note; sprintf( "fminopt.Algorithm:\t%s", string( oc.fminopt.Algorithm ) ) ];
    note = [ note; sprintf( "fminopt.FiniteDifferenceStepSize:\t%s", string( oc.fminopt.FiniteDifferenceStepSize ) ) ];
    note = [ note; sprintf( "fminopt.FunctionTolerance:\t%s", string( oc.fminopt.FunctionTolerance ) ) ];
    note = [ note; sprintf( "fminopt.ConstraintTolerance:\t%s", string( oc.fminopt.ConstraintTolerance ) ) ];
    if isfield( oc.fminopt, 'OutputFcn' )
        note = [ note; sprintf( "fminopt.OutputFcn:\t%s", func2str( oc.fminopt.OutputFcn ) ) ];
    end
    if matches( oc.fminopt.Algorithm, "interior-point", 'ignorecase', true )
        note = [ note; sprintf( "fminopt.BarrierParamUpdate:\t%s", string( oc.fminopt.BarrierParamUpdate ) ) ];
        note = [ note; sprintf( "fminopt.EnableFeasibilityMode:\t%s", string( oc.fminopt.EnableFeasibilityMode ) ) ];
        note = [ note; sprintf( "fminopt.HessianApproximation:\t%s", string( oc.fminopt.HessianApproximation ) ) ];
        note = [ note; sprintf( "fminopt.HonorBounds:\t%s", string( oc.fminopt.HonorBounds ) ) ];
        note = [ note; sprintf( "fminopt.InitBarrierParam:\t%s", string( oc.fminopt.InitBarrierParam ) ) ];
        note = [ note; sprintf( "fminopt.InitTrustRegionRadius:\t%s", string( oc.fminopt.InitTrustRegionRadius ) ) ];
        note = [ note; sprintf( "fminopt.MaxProjCGIter:\t%s", string( oc.fminopt.MaxProjCGIter ) ) ];
        note = [ note; sprintf( "fminopt.ObjectiveLimit:\t%s", string( oc.fminopt.ObjectiveLimit ) ) ];
        note = [ note; sprintf( "fminopt.ScaleProblem:\t%s", string( oc.fminopt.ScaleProblem ) ) ];
        note = [ note; sprintf( "fminopt.SubproblemAlgorithm:\t%s", string( oc.fminopt.SubproblemAlgorithm ) ) ];
        note = [ note; sprintf( "fminopt.TolProjCG:\t%s", string( oc.fminopt.TolProjCG ) ) ];
        note = [ note; sprintf( "fminopt.TolProjCGAbs:\t%s", string( oc.fminopt.TolProjCGAbs ) ) ];
    elseif matches( oc.fminopt.Algorithm, "active-set", 'ignorecase', true )
        note = [ note; sprintf( "fminopt.MaxSQPIter:\t%s", string( oc.fminopt.MaxSQPIter ) ) ];
        note = [ note; sprintf( "fminopt.RelLineSrchBnd:\t%s", string( oc.fminopt.RelLineSrchBnd ) ) ];
        note = [ note; sprintf( "fminopt.RelLineSrchBndDuration:\t%s", string( oc.fminopt.RelLineSrchBndDuration ) ) ];
        note = [ note; sprintf( "fminopt.TolConSQP:\t%s", string( oc.fminopt.TolConSQP ) ) ];
    elseif matches( oc.fminopt.Algorithm, "sqp", 'ignorecase', true )
        note = [ note; sprintf( "fminopt.ObjectiveLimit:\t%s", string( oc.fminopt.ObjectiveLimit ) ) ];
        note = [ note; sprintf( "fminopt.ScaleProblem:\t%s", string( oc.fminopt.ScaleProblem ) ) ];
    end
end
if matches( oc.optType, "ga", 'ignorecase', true )
    note = [ note; "" ];
    note = [ note; "ga options:" ];
    note = [ note; sprintf( "gaopt.UseParallel:\t%s", string( oc.gaopt.UseParallel ) ) ];
    note = [ note; sprintf( "gaopt.MaxGenerations:\t%s", string( oc.gaopt.MaxGenerations ) ) ];
    note = [ note; sprintf( "gaopt.MaxStallGenerations:\t%s", string( oc.gaopt.MaxStallGenerations ) ) ];
    note = [ note; sprintf( "gaopt.PopulationSize:\t%s", string( oc.gaopt.PopulationSize ) ) ];
    note = [ note; sprintf( "gaopt.MaxTime:\t%s", string( oc.gaopt.MaxTime ) ) ];
    note = [ note; sprintf( "gaopt.MaxStallTime:\t%s", string( oc.gaopt.MaxStallTime ) ) ];
    note = [ note; sprintf( "gaopt.CrossoverFraction:\t%s", string( oc.gaopt.CrossoverFraction ) ) ];
    note = [ note; sprintf( "gaopt.MutationFcn:\t%s", string( oc.gaopt.MutationFcn ) ) ];
    note = [ note; sprintf( "gaopt.CreationFcn:\t%s", string( oc.gaopt.CreationFcn ) ) ];
    note = [ note; sprintf( "gaopt.CrossoverFcn:\t%s", string( oc.gaopt.CrossoverFcn ) ) ];
    note = [ note; sprintf( "gaopt.SelectionFcn:\t%s", string( oc.gaopt.SelectionFcn ) ) ];
    note = [ note; ...
        strcat( sprintf( "gaopt.InitialPopulationRange:\t[" ), sprintf( " %g ", oc.gaopt.InitialPopulationRange ), sprintf("]") ) ];
    if isfield( oc, 'rngstate')
        note = [ note; sprintf( "gaopt.Type:\t%s", string( oc.gaopt.Type ) ) ];
        note = [ note; sprintf( "gaopt.Seed:\t%s", string( oc.gaopt.Seed ) ) ];
        note = [ note; sprintf( "gaopt.State:\t%s", string( oc.gaopt.State ) ) ];
    end
end
if matches( oc.optType, "ipopt", 'ignorecase', true )
    note = [ note; "" ];
    note = [ note; "ipopt options:" ];
    note = [ note; sprintf( "ipopt.options.ipopt.print_level:\t%s", string( oc.ipopt.options.ipopt.print_level ) ) ];
    if isfield( oc.ipopt.options.ipopt, 'tol' )
        note = [ note; sprintf( "ipopt.options.ipopt.tol:\t%s", string( oc.ipopt.options.ipopt.tol ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'constr_viol_tol' )
        note = [ note; sprintf( "ipopt.options.ipopt.constr_viol_tol:\t%s", string( oc.ipopt.options.ipopt.constr_viol_tol ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'dual_inf_tol' )
        note = [ note; sprintf( "ipopt.options.ipopt.dual_inf_tol:\t%s", string( oc.ipopt.options.ipopt.dual_inf_tol ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'compl_inf_tol' )
        note = [ note; sprintf( "ipopt.options.ipopt.compl_inf_tol:\t%s", string( oc.ipopt.options.ipopt.compl_inf_tol ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'print_user_options' )
        note = [ note; sprintf( "ipopt.options.ipopt.print_user_options:\t%s", string( oc.ipopt.options.ipopt.print_user_options ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'max_iter' )
        note = [ note; sprintf( "ipopt.options.ipopt.max_iter:\t%s", string( oc.ipopt.options.ipopt.max_iter ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'linear_solver' )
        note = [ note; sprintf( "ipopt.options.ipopt.linear_solver:\t%s", string( oc.ipopt.options.ipopt.linear_solver ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'mu_strategy' )
        note = [ note; sprintf( "ipopt.options.ipopt.mu_strategy:\t%s", string( oc.ipopt.options.ipopt.mu_strategy ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'mu_init' )
        note = [ note; sprintf( "ipopt.options.ipopt.mu_init:\t%s", string( oc.ipopt.options.ipopt.mu_init ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'mu_min' )
        note = [ note; sprintf( "ipopt.options.ipopt.mu_min:\t%s", string( oc.ipopt.options.ipopt.mu_min ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'mu_max' )
        note = [ note; sprintf( "ipopt.options.ipopt.mu_max:\t%s", string( oc.ipopt.options.ipopt.mu_max ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'adaptive_mu_globalization' )
        note = [ note; sprintf( "ipopt.options.ipopt.adaptive_mu_globalization:\t%s", string( oc.ipopt.options.ipopt.adaptive_mu_globalization ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'barrier_tol_factor' )
        note = [ note; sprintf( "ipopt.options.ipopt.barrier_tol_factor:\t%s", string( oc.ipopt.options.ipopt.barrier_tol_factor ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'alpha_for_y' )
        note = [ note; sprintf( "ipopt.options.ipopt.alpha_for_y:\t%s", string( oc.ipopt.options.ipopt.alpha_for_y ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'recalc_y' )
        note = [ note; sprintf( "ipopt.options.ipopt.recalc_y:\t%s", string( oc.ipopt.options.ipopt.recalc_y ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'recalc_y_feas_tol' )
        note = [ note; sprintf( "ipopt.options.ipopt.recalc_y_feas_tol:\t%s", string( oc.ipopt.options.ipopt.recalc_y_feas_tol ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'max_soc' )
        note = [ note; sprintf( "ipopt.options.ipopt.max_soc:\t%s", string( oc.ipopt.options.ipopt.max_soc ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'alpha_min_frac' )
        note = [ note; sprintf( "ipopt.options.ipopt.alpha_min_frac:\t%s", string( oc.ipopt.options.ipopt.alpha_min_frac ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'alpha_red_factor' )
        note = [ note; sprintf( "ipopt.options.ipopt.alpha_red_factor:\t%s", string( oc.ipopt.options.ipopt.alpha_red_factor ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'kappa_sigma' )
        note = [ note; sprintf( "ipopt.options.ipopt.kappa_sigma:\t%s", string( oc.ipopt.options.ipopt.kappa_sigma ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'limited_memory_max_history' )
        note = [ note; sprintf( "ipopt.options.ipopt.limited_memory_max_history:\t%s", string( oc.ipopt.options.ipopt.limited_memory_max_history ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'nlp_scaling_max_gradient' )
        note = [ note; sprintf( "ipopt.options.ipopt.nlp_scaling_max_gradient:\t%s", string( oc.ipopt.options.ipopt.nlp_scaling_max_gradient ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'mu_linear_decrease_factor' )
        note = [ note; sprintf( "ipopt.options.ipopt.mu_linear_decrease_factor:\t%s", string( oc.ipopt.options.ipopt.mu_linear_decrease_factor  ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'mu_superlinear_decrease_power' )
        note = [ note; sprintf( "ipopt.options.ipopt.mu_superlinear_decrease_power:\t%s", string( oc.ipopt.options.ipopt.mu_superlinear_decrease_power  ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'max_wall_time' )
        note = [ note; sprintf( "ipopt.options.ipopt.max_wall_time:\t%s", string( oc.ipopt.options.ipopt.max_wall_time  ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'bound_push' )
        note = [ note; sprintf( "ipopt.options.ipopt.bound_push:\t%s", string( oc.ipopt.options.ipopt.bound_push  ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'bound_frac' )
        note = [ note; sprintf( "ipopt.options.ipopt.bound_frac:\t%s", string( oc.ipopt.options.ipopt.bound_frac  ) ) ];
    end
    if isfield( oc.ipopt.options.ipopt, 'bound_relax_factor' )
        note = [ note; sprintf( "ipopt.options.ipopt.bound_relax_factor:\t%s", string( oc.ipopt.options.ipopt.bound_relax_factor  ) ) ];
    end
end
if matches( oc.optType, "snopt", 'ignorecase', true )
    note = [ note; "" ];
    note = [ note; "snopt options:" ];
    note = [ note; sprintf( "snopt.start:\t%s", string( oc.snopt.start ) ) ];
    note = [ note; sprintf( "snopt.majorfeasibilitytolerance:\t%s", string( oc.snopt.majorfeasibilitytolerance ) ) ];
    note = [ note; sprintf( "snopt.majoroptimalitytolerance:\t%s", string( oc.snopt.majoroptimalitytolerance ) ) ];
    note = [ note; sprintf( "snopt.minorfeasibilitytolerance:\t%s", string( oc.snopt.minorfeasibilitytolerance ) ) ];
    note = [ note; sprintf( "snopt.minoroptimalitytolerance:\t%s", string( oc.snopt.minoroptimalitytolerance ) ) ];
    note = [ note; sprintf( "snopt.scaleoption:\t%s", string( oc.snopt.scaleoption ) ) ];
    note = [ note; sprintf( "snopt.minoriterationslimit:\t%s", string( oc.snopt.minoriterationslimit ) ) ];
    note = [ note; sprintf( "snopt.majorsteplimit:\t%s", string( oc.snopt.majorsteplimit ) ) ];
    note = [ note; sprintf( "snopt.derivativelevel:\t%s", string( oc.snopt.derivativelevel ) ) ];
    note = [ note; sprintf( "snopt.majoriterationslimit:\t%s", string( oc.snopt.majoriterationslimit ) ) ];
    if isfield( oc.snopt, 'timelimit' ) 
        note = [ note; sprintf( "snopt.timelimit:\t%s", string( oc.snopt.timelimit ) ) ];
    end
end
if matches( oc.optType, "nlopt", 'ignorecase', true )
    note = [ note; "" ];
    note = [ note; "nlopt options:" ];
    note = [ note; sprintf( "nlopt.algorithm:\t%s", string( oc.nlopt.algorithm ) ) ];
    note = [ note; sprintf( "nlopt.ConstraintTolerance:\t%s", string( oc.nlopt.ConstraintTolerance ) ) ];
    note = [ note; sprintf( "nlopt.stopval:\t%s", string( oc.nlopt.stopval ) ) ];
    note = [ note; sprintf( "nlopt.ftol_rel:\t%s", string( oc.nlopt.ftol_rel ) ) ];
    note = [ note; sprintf( "nlopt.ftol_abs:\t%s", string( oc.nlopt.ftol_abs ) ) ];
    note = [ note; sprintf( "nlopt.xtol_rel:\t%s", string( oc.nlopt.xtol_rel ) ) ];
    note = [ note; sprintf( "nlopt.maxtime:\t%s", string( oc.nlopt.maxtime ) ) ];
    note = [ note; sprintf( "nlopt.maxeval:\t%s", string( oc.nlopt.maxeval ) ) ];
end
note = [ note; "" ];
note = [ note; "" ];

%% Output Metrics
note = [ note; "Output Metrics:" ];

note = [ note; "" ];
note = [ note; "opt metrics:" ];
if isfield( opt, 'initOptTime' )
    note = [ note; sprintf( "opt.initOptTime:\t%g", opt.initOptTime ) ];
end
note = [ note; sprintf( "opt.output.fval:\t%g", opt.output.fval ) ];
note = [ note; sprintf( "opt.output.exitflag:\t%i", opt.output.exitflag ) ];
if isfield( opt.output, 'output' )
    if isfield( opt.output.output, 'iterations' )
        note = [ note; sprintf( "opt.output.iterations:\t%g", opt.output.output.iterations ) ];
    end
    if isfield( opt.output.output, 'funcCount' )
        note = [ note; sprintf( "opt.output.funcCount:\t%g", opt.output.output.funcCount ) ];
    end
end
note = [ note; sprintf( "opt.optTime:\t%g", opt.optTime ) ];
note = [ note; sprintf( "opt.numPos:\t%i", opt.numPos ) ];
note = [ note; sprintf( "opt.numVars:\t%i", opt.numVars ) ];
note = [ note; sprintf( "opt.useGPU:\t%s", string( opt.useGPU ) ) ];

note = [ note; "" ];
note = [ note; "opt constraint metrics:" ];
constvalsList = keys( opt.constraintvals );
for kk = 1:length( constvalsList )
    note = [ note; ...
        sprintf( "%s:\t%06.3f%%", constvalsList( kk ), opt.constraintvals( constvalsList( kk ) ) * 100 ) ]; %#ok
end

note = [ note; "" ];
note = [ note; "waveform metrics:" ];
note = [ note; sprintf( "breal_phasor_slew_max:\t%g", valtrain.breal_phasor_slew_max ) ];
note = [ note; sprintf( "breal_phasor_accel_max:\t%g", valtrain.breal_phasor_accel_max ) ];
note = [ note; sprintf( "bimag_phasor_slew_max:\t%g", valtrain.bimag_phasor_slew_max ) ];
note = [ note; sprintf( "bimag_phasor_accel_max:\t%g", valtrain.bimag_phasor_accel_max ) ];
note = [ note; sprintf( "Grad_slew_max:\t%g", valtrain.Grad_slew_max ) ];
note = [ note; sprintf( "Grad_accel_max:\t%g", valtrain.Grad_accel_max ) ];
note = [ note; sprintf( "Shim_slew_max:\t%g", valtrain.Shim_slew_max ) ];
note = [ note; sprintf( "Shim_accel_max:\t%g", valtrain.Shim_accel_max ) ];

if pulse.sliceSelective

    note = [ note; "" ];
    note = [ note; "valtrain metrics:" ];
    note = [ note; sprintf( "valtrain.inSliceRipple:\t%g", mean( valtrain.inSliceRipple, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.outOfSliceRipple:\t%g", mean( valtrain.outOfSliceRipple, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.inSliceRipplePercent:\t%g", mean( valtrain.inSliceRipplePercent, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.outOfSliceRipplePercent:\t%g", mean( valtrain.outOfSliceRipplePercent, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.transitionWidth:\t%g", mean( valtrain.transitionWidth, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.FWHM:\t%g", mean( valtrain.FWHM, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.magNRMSEThruPlane_FA:\t%g", mean( valtrain.magNRMSEThruPlane_FA, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.magNRMSEThruPlane_compFA:\t%g", mean( valtrain.magNRMSEThruPlane_compFA, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.magNRMSEInPlane_Mxy:\t%g", mean( valtrain.magNRMSEInPlane_Mxy, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.magNRMSEInPlane_FA:\t%g", mean( valtrain.magNRMSEInPlane_FA, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.meanFAInPlane:\t%g", mean( valtrain.meanFAInPlane, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.stdFAInPlane:\t%g", mean( valtrain.stdFAInPlane, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.nMaxMinInPlane:\t%g", mean( valtrain.nMaxMinInPlane, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.magCoeffVarInPlane:\t%g", mean( valtrain.magCoeffVarInPlane, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.magNRMSEInPlaneCentralSlice_Mxy:\t%g", mean( valtrain.magNRMSEInPlaneCentralSlice_Mxy, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.NRMSEInPlaneCentralSlice_Mxy:\t%g", mean( valtrain.NRMSEInPlaneCentralSlice_Mxy, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.magNRMSEInPlaneCentralSlice_FA:\t%g", mean( valtrain.magNRMSEInPlaneCentralSlice_FA, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.meanFAInPlaneCentralSlice:\t%g", mean( valtrain.meanFAInPlaneCentralSlice, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.stdFAInPlaneCentralSlice:\t%g", mean( valtrain.stdFAInPlaneCentralSlice, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.nMaxMinInPlaneCentralSlice:\t%g", mean( valtrain.nMaxMinInPlaneCentralSlice, 'all' ) ) ];
    note = [ note; sprintf( "valtrain.magCoeffVarInPlaneCentralSlice:\t%g", mean( valtrain.magCoeffVarInPlaneCentralSlice, 'all' ) ) ];

else
    note = [ note; "" ];
    note = [ note; "valtrain metrics:" ];
    note = [ note; sprintf( "valtrain.magNRMSE_Mxy:\t%g", mean( valtrain.magNRMSE_Mxy ) ) ];
    note = [ note; sprintf( "valtrain.NRMSE_Mxy:\t%g", mean( valtrain.NRMSE_Mxy ) ) ];
    note = [ note; sprintf( "valtrain.magNRMSE_FA:\t%g", mean( valtrain.magNRMSE_FA ) ) ];
    note = [ note; sprintf( "valtrain.meanFA:\t%g", mean( valtrain.meanFA ) ) ];
    note = [ note; sprintf( "valtrain.stdFA:\t%g", mean( valtrain.stdFA ) ) ];
    note = [ note; sprintf( "valtrain.nMaxMin:\t%g", mean( valtrain.nMaxMin ) ) ];
    note = [ note; sprintf( "valtrain.magCoeffVar:\t%g", mean( valtrain.magCoeffVar ) ) ];

    if ( valtrain.numXYCoils == 1 ) && valtrain.performBCHPcomp
        note = [ note; sprintf( "valtrain.BCHP_magNRMSE_Mxy:\t%g", mean( valtrain.BCHP_magNRMSE_Mxy ) ) ];
        note = [ note; sprintf( "valtrain.BCHP_NRMSE_Mxy:\t%g", mean( valtrain.BCHP_NRMSE_Mxy ) ) ];
        note = [ note; sprintf( "valtrain.BCHP_magNRMSE_FA:\t%g", mean( valtrain.BCHP_magNRMSE_FA ) ) ];
        note = [ note; sprintf( "valtrain.BCHP_meanFA:\t%g", mean( valtrain.BCHP_meanFA ) ) ];
        note = [ note; sprintf( "valtrain.BCHP_stdFA:\t%g", mean( valtrain.BCHP_stdFA ) ) ];
        note = [ note; sprintf( "valtrain.BCHP_nMaxMin:\t%g", mean( valtrain.BCHP_nMaxMin ) ) ];
        note = [ note; sprintf( "valtrain.BCHP_magCoeffVar:\t%g", mean( valtrain.BCHP_magCoeffVar ) ) ];
        note = [ note; ...
            strcat( sprintf( "valtrain.BCHP_BCmag:\t[" ), sprintf( " %g ", valtrain.BCHP_BCmag ), sprintf("]") ) ];
        note = [ note; ...
            strcat( sprintf( "valtrain.BCHP_RFPower:\t[" ), sprintf( " %g ", valtrain.BCHP_RFPower ), sprintf("]") ) ];
    end
end

note = [ note; sprintf( "valtrain.simTime:\t%g", valtrain.simTime ) ];
note = [ note; sprintf( "valtrain.totalRFPower:\t%g", valtrain.totalRFPower ) ];
note = [ note; sprintf( "valtrain.maxRFPower:\t%g", valtrain.maxRFPower ) ];

if isfield( valtrain, 'VOPs' )
    note = [ note; sprintf( "valtrain.peakLocalSAR:\t%g", valtrain.peakLocalSAR ) ];
    note = [ note; sprintf( "valtrain.peakAvgLocalSAR:\t%g", valtrain.peakAvgLocalSAR ) ];
end
if isfield( valtrain, 'QGlobal' )
    note = [ note; sprintf( "valtrain.peakGlobalSAR:\t%g", valtrain.peakGlobalSAR ) ];
    note = [ note; sprintf( "valtrain.avgGlobalSAR:\t%g", valtrain.avgGlobalSAR ) ];
end

if ~isempty( valtest )
    if pulse.sliceSelective

        note = [ note; "" ];
        note = [ note; "valtest metrics:" ];
        note = [ note; sprintf( "valtest.inSliceRipple:\t%g", mean( valtest.inSliceRipple, 'all' ) ) ];
        note = [ note; sprintf( "valtest.outOfSliceRipple:\t%g", mean( valtest.outOfSliceRipple, 'all' ) ) ];
        note = [ note; sprintf( "valtest.inSliceRipplePercent:\t%g", mean( valtest.inSliceRipplePercent, 'all' ) ) ];
        note = [ note; sprintf( "valtest.outOfSliceRipplePercent:\t%g", mean( valtest.outOfSliceRipplePercent, 'all' ) ) ];
        note = [ note; sprintf( "valtest.transitionWidth:\t%g", mean( valtest.transitionWidth, 'all' ) ) ];
        note = [ note; sprintf( "valtest.FWHM:\t%g", mean( valtest.FWHM, 'all' ) ) ];
        note = [ note; sprintf( "valtest.magNRMSEThruPlane_FA:\t%g", mean( valtest.magNRMSEThruPlane_FA, 'all' ) ) ];
        note = [ note; sprintf( "valtest.magNRMSEThruPlane_compFA:\t%g", mean( valtest.magNRMSEThruPlane_compFA, 'all' ) ) ];
        note = [ note; sprintf( "valtest.magNRMSEInPlane_Mxy:\t%g", mean( valtest.magNRMSEInPlane_Mxy, 'all' ) ) ];
        note = [ note; sprintf( "valtest.magNRMSEInPlane_FA:\t%g", mean( valtest.magNRMSEInPlane_FA, 'all' ) ) ];
        note = [ note; sprintf( "valtest.meanFAInPlane:\t%g", mean( valtest.meanFAInPlane, 'all' ) ) ];
        note = [ note; sprintf( "valtest.stdFAInPlane:\t%g", mean( valtest.stdFAInPlane, 'all' ) ) ];
        note = [ note; sprintf( "valtest.nMaxMinInPlane:\t%g", mean( valtest.nMaxMinInPlane, 'all' ) ) ];
        note = [ note; sprintf( "valtest.magCoeffVarInPlane:\t%g", mean( valtest.magCoeffVarInPlane, 'all' ) ) ];
        note = [ note; sprintf( "valtest.magNRMSEInPlaneCentralSlice_Mxy:\t%g", mean( valtest.magNRMSEInPlaneCentralSlice_Mxy, 'all' ) ) ];
        note = [ note; sprintf( "valtest.NRMSEInPlaneCentralSlice_Mxy:\t%g", mean( valtest.NRMSEInPlaneCentralSlice_Mxy, 'all' ) ) ];
        note = [ note; sprintf( "valtest.magNRMSEInPlaneCentralSlice_FA:\t%g", mean( valtest.magNRMSEInPlaneCentralSlice_FA, 'all' ) ) ];
        note = [ note; sprintf( "valtest.meanFAInPlaneCentralSlice:\t%g", mean( valtest.meanFAInPlaneCentralSlice, 'all' ) ) ];
        note = [ note; sprintf( "valtest.stdFAInPlaneCentralSlice:\t%g", mean( valtest.stdFAInPlaneCentralSlice, 'all' ) ) ];
        note = [ note; sprintf( "valtest.nMaxMinInPlaneCentralSlice:\t%g", mean( valtest.nMaxMinInPlaneCentralSlice, 'all' ) ) ];
        note = [ note; sprintf( "valtest.magCoeffVarInPlaneCentralSlice:\t%g", mean( valtest.magCoeffVarInPlaneCentralSlice, 'all' ) ) ];

    else
        note = [ note; "" ];
        note = [ note; "valtest metrics:" ];
        note = [ note; sprintf( "valtest.magNRMSE_Mxy:\t%g", mean( valtest.magNRMSE_Mxy ) ) ];
        note = [ note; sprintf( "valtest.NRMSE_Mxy:\t%g", mean( valtest.NRMSE_Mxy ) ) ];
        note = [ note; sprintf( "valtest.magNRMSE_FA:\t%g", mean( valtest.magNRMSE_FA ) ) ];
        note = [ note; sprintf( "valtest.meanFA:\t%g", mean( valtest.meanFA ) ) ];
        note = [ note; sprintf( "valtest.stdFA:\t%g", mean( valtest.stdFA ) ) ];
        note = [ note; sprintf( "valtest.nMaxMin:\t%g", mean( valtest.nMaxMin ) ) ];
        note = [ note; sprintf( "valtest.magCoeffVar:\t%g", mean( valtest.magCoeffVar ) ) ];

        if ( valtest.numXYCoils == 1 ) && valtest.performBCHPcomp
            note = [ note; sprintf( "valtest.BCHP_magNRMSE_Mxy:\t%g", mean( valtest.BCHP_magNRMSE_Mxy ) ) ];
            note = [ note; sprintf( "valtest.BCHP_NRMSE_Mxy:\t%g", mean( valtest.BCHP_NRMSE_Mxy ) ) ];
            note = [ note; sprintf( "valtest.BCHP_magNRMSE_FA:\t%g", mean( valtest.BCHP_magNRMSE_FA ) ) ];
            note = [ note; sprintf( "valtest.BCHP_meanFA:\t%g", mean( valtest.BCHP_meanFA ) ) ];
            note = [ note; sprintf( "valtest.BCHP_stdFA:\t%g", mean( valtest.BCHP_stdFA ) ) ];
            note = [ note; sprintf( "valtest.BCHP_nMaxMin:\t%g", mean( valtest.BCHP_nMaxMin ) ) ];
            note = [ note; sprintf( "valtest.BCHP_magCoeffVar:\t%g", mean( valtest.BCHP_magCoeffVar ) ) ];
            note = [ note; ...
                strcat( sprintf( "valtest.BCHP_BCmag:\t[" ), sprintf( " %g ", valtest.BCHP_BCmag ), sprintf("]") ) ];
            note = [ note; ...
                strcat( sprintf( "valtest.BCHP_RFPower:\t[" ), sprintf( " %g ", valtest.BCHP_RFPower ), sprintf("]") ) ];
        end
    end

    note = [ note; sprintf( "valtest.simTime:\t%g", valtest.simTime ) ];
    note = [ note; sprintf( "valtest.totalRFPower:\t%g", valtest.totalRFPower ) ];
    note = [ note; sprintf( "valtest.maxRFPower:\t%g", valtest.maxRFPower ) ];

    if isfield( valtest, 'VOPs' )
        note = [ note; sprintf( "valtest.peakLocalSAR:\t%g", valtest.peakLocalSAR ) ];
        note = [ note; sprintf( "valtest.peakAvgLocalSAR:\t%g", valtest.peakAvgLocalSAR ) ];
    end
    if isfield( valtest, 'QGlobal' )
        note = [ note; sprintf( "valtest.peakGlobalSAR:\t%g", valtest.peakGlobalSAR ) ];
        note = [ note; sprintf( "valtest.avgGlobalSAR:\t%g", valtest.avgGlobalSAR ) ];
    end

end

note = [ note; "" ];
note = [ note; "" ];

end