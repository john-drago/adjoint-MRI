function oc = getOC3D( ocType, st )

%% OC preliminaries
oc = struct; % define opt control struct
oc.binVis = st.binVis; % whether or not to make images visible
oc.saveResult = st.saveResult; % whether or not to save optimization results
oc.saveDir = st.saveDir; % location of where to save optimization results
if isfield( st, 'note' )
    oc.note = st.note;
end
oc.saveFileRecord = st.saveFileRecord;
oc.currFile = st.currFile;
oc.useGPU = st.useGPU;
oc.trackConvergence = st.trackConvergence;
oc.trackDecisionVariables = st.trackDecisionVariables;
if isfield( st, 'spacingTrackDecisionVariables' )
    oc.spacingTrackDecisionVariables = st.spacingTrackDecisionVariables;
end
oc.ensureFeasibleStart = st.ensureFeasibleStart;

%% OC spacing

oc.opt_di = st.opt_di .* ones( 3, 1 ); % dx, dy, dz for optimization
oc.opt_dt = st.opt_dt; % time in seconds

oc.val_di = st.val_di .* ones( 3, 1 ); % dx, dy, dz for validation
oc.val_dt = st.val_dt; % time in seconds

% optimization specific parameters
if strcmpi( ocType, "active-set" ) || strcmpi( ocType, "interior-point" )
    oc.optType = "fmin";

    optDisplay = "iter";

    oc.fminopt = optimoptions( "fmincon" );
    oc.fminopt.Display = optDisplay;
    oc.fminopt.UseParallel = st.useParallel;
    oc.fminopt.MaxFunctionEvaluations = inf;
    oc.fminopt.FiniteDifferenceStepSize = 1e-7;
    oc.fminopt.FunctionTolerance = 1e-6; % doesn't matter for SQP
    oc.fminopt.ConstraintTolerance = 1e-6;
    oc.fminopt.SpecifyObjectiveGradient = true;
    oc.fminopt.SpecifyConstraintGradient = true;
    oc.fminopt.CheckGradients = false;

    oc.fminopt.MaxIterations = 300;

end

if strcmpi( ocType, "active-set" )

    oc.fminopt.Algorithm = 'active-set';
    oc.fminopt.StepTolerance = 1e-8;

    oc.fminopt.RelLineSrchBnd = 5e-3;
    oc.fminopt.RelLineSrchBndDuration = inf;
    % oc.fminopt.TolConSQP = 1e-8;
    
    oc.fminopt.MaxIterations = 300;
    % oc.fminopt.MaxIterations = 3;

elseif strcmpi( ocType, "sqp" )
    oc.fminopt.Algorithm = 'sqp';
    oc.fminopt.StepTolerance = 1e-6;

    oc.fminopt.MaxIterations = 300;
    % oc.fminopt.MaxIterations = 3;

elseif strcmpi( ocType, "interior-point" )

    oc.fminopt.Algorithm = 'interior-point';
    oc.fminopt.BarrierParamUpdate = 'predictor-corrector';
    oc.fminopt.StepTolerance = 1e-10;
    oc.fminopt.HessianApproximation = 'bfgs';
    oc.fminopt.InitBarrierParam = 1e-5;
    % oc.fminopt.EnableFeasibilityMode = true;
    oc.fminopt.ScaleProblem = true;
    oc.fminopt.SubproblemAlgorithm = 'cg';
    oc.fminopt.TolProjCG = 1e-8;
    oc.fminopt.TolProjCGAbs = 1e-8;

    oc.fminopt.MaxIterations = 300;
    % oc.fminopt.MaxIterations = 3;

elseif strcmpi( ocType, "ipopt" )

    oc.optType = 'ipopt';
    oc.ipopt = struct;
    oc.ipopt.options = struct;
    oc.ipopt.options.ipopt = struct;
    oc.ipopt.options.ipopt.print_user_options = 'yes';
    % oc.ipopt.options.ipopt.print_options_documentation = 'yes';
    % oc.ipopt.options.ipopt.print_advanced_options = 'yes';
    oc.ipopt.options.ipopt.print_level = 5;

    oc.ipopt.options.ipopt.mu_strategy = 'adaptive';
    oc.ipopt.options.ipopt.mu_min = 1e-11;
    oc.ipopt.options.ipopt.mu_max = 1e-5;
    oc.ipopt.options.ipopt.adaptive_mu_globalization = 'kkt-error';

    oc.ipopt.options.ipopt.tol = 1e-6;
    oc.ipopt.options.ipopt.constr_viol_tol = 1e-6;
    oc.ipopt.options.ipopt.dual_inf_tol = 1e-6;
    oc.ipopt.options.ipopt.compl_inf_tol = 1e-6;

    bound_push_frac = 1e-4;
    oc.ipopt.options.ipopt.bound_push = bound_push_frac;
    oc.ipopt.options.ipopt.bound_frac = bound_push_frac;
    oc.ipopt.options.ipopt.bound_relax_factor = 1e-4;

    oc.ipopt.options.ipopt.alpha_for_y = 'safer-min-dual-infeas';
    oc.ipopt.options.ipopt.alpha_min_frac = 0.01;
    oc.ipopt.options.ipopt.max_soc = 10;
    oc.ipopt.options.ipopt.recalc_y = 'yes';
    oc.ipopt.options.ipopt.recalc_y_feas_tol = 1e-2;

    % oc.ipopt.options.ipopt.kappa_sigma = 1e+6;

    oc.ipopt.options.ipopt.nlp_scaling_method = 'gradient-based';
    oc.ipopt.options.ipopt.nlp_scaling_max_gradient = 100;
    oc.ipopt.options.ipopt.nlp_scaling_min_value = 1e-8;
    oc.ipopt.options.ipopt.linear_solver = 'mumps';

    oc.ipopt.options.ipopt.max_iter = 300;
    % oc.ipopt.options.ipopt.max_iter = 3;

end

%% Final process and output
oc = determineParallelWorkers( st.numWorkers, oc );

oc %#ok print oc information

end