function [ oc, pulse ] = getOC2D( ocType, pulse, fields, st )

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
oc.opt_di = struct;

% withinSlice_di = 0.001;
% outsideSlice_di = 0.015;
% intraSlice_di = 0.025;
% outsideSliceExt = 1.0;
% % outsideSlice_di_start = 0.002;
% % outsideSlice_di_end = 0.005;

withinSlice_di = 0.00025;
outsideSlice_di = 0.002;
intraSlice_di = 0.010;
outsideSliceExt = 1.50;
% outsideSlice_di_start = 0.002;
% outsideSlice_di_end = 0.005;

oc.opt_di.withinSlice_di = withinSlice_di;
oc.opt_di.outsideSlice_di = outsideSlice_di;
oc.opt_di.intraSlice_di = intraSlice_di;
oc.opt_di.outsideSliceExt = outsideSliceExt;
% oc.opt_di.outsideSlice_di_start = outsideSlice_di_start;
% oc.opt_di.outsideSlice_di_end = outsideSlice_di_end;

pulse.sliceBounds = determineSliceBounds( pulse.sliceLocation, pulse.sliceThickness );

[ oc.opt_di.z, pulse.sliceBounds ] = getOptimizationDiscretizationLinspaceSliceSelect(...
    withinSlice_di, outsideSlice_di, outsideSliceExt,...
    pulse.sliceLocation, pulse.sliceThickness, pulse.sliceBounds, [ fields.z(1), fields.z(end) ] );

% [ oc.opt_di.z, pulse.sliceBounds ] = getOptimizationDiscretizationLogspaceSliceSelect(...
%     withinSlice_di, outsideSlice_di_start, outsideSlice_di_end,...
%     outsideSliceExt, pulse.sliceLocation, pulse.sliceThickness, pulse.sliceBounds, [ fields.z(1), fields.z(end) ] );

val_di = st.val_di;

pulse.targSt.x = fields.x(1) : val_di : fields.x(end);
pulse.targSt.y = fields.y(1) : val_di : fields.y(end);
pulse.targSt.z = unique( [...
    fields.z(1) : val_di : fields.z(end),...
    oc.opt_di.z ] );

oc.opt_di.x = fields.x(1) : intraSlice_di : fields.x(end);
oc.opt_di.y = fields.y(1) : intraSlice_di : fields.y(end);

% x_const = [];
% y_const = [];
% z_const = unique( [...
%     pulse.sliceBounds( : );...
%     pulse.sliceBounds( :, 1 ) - withinSlice_di;...
%     pulse.sliceBounds( :, 2 ) + withinSlice_di  ] );
% [ oc.opt_di.x, oc.opt_di.y, oc.opt_di.z ] =  perturbOptimizationPoints(...
%     oc.opt_di.x, oc.opt_di.y, oc.opt_di.z, x_const, y_const, z_const );

oc.opt_dt = pulse.dt; % time in seconds

pulse.targSt.MtargVals = [ sind( pulse.targFAVal ); 0; cosd( pulse.targFAVal ) ];
pulse.targSt = generateMTargSliceSelect( pulse.targSt, pulse );

oc.val_di = zeros( 3, 1 );
oc.val_di( 1:3, 1 ) = val_di;
oc.val_dt = pulse.dt; % time in seconds

%% OC solver type
% optimization specific parameters
if matches( ocType, [ "active-set", "interior-point", "sqp" ], 'ignorecase', true )

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

end

if strcmpi( ocType, "active-set" )

    oc.fminopt.Algorithm = 'active-set';
    oc.fminopt.StepTolerance = 1e-8;
    oc.fminopt.RelLineSrchBnd = 1.50e-4;
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

    oc.ipopt.options.ipopt.max_iter = 300;

end

%% Final process and output
oc = determineParallelWorkers( st.numWorkers, oc );

oc %#ok print oc information

end