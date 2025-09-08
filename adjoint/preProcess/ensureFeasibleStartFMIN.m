function [ pSc0Opt, feasStartOutput ] = ensureFeasibleStartFMIN(...
    pSc0, opt, efsSettings )

if nargin < 3
    efsSettings = struct;
end

if ~isfield( efsSettings, 'relLineSrchBnd' )
    efsSettings.relLineSrchBnd = 5e-5;
end
if ~isfield( efsSettings, 'maxTime' )
    efsSettings.maxTime = 63 * 60; % seconds
end
if ~isfield( efsSettings, 'scaleAtEnd' )
    efsSettings.scaleAtEnd = false;
end
if ~isfield( efsSettings, 'usenonlcon' )
    efsSettings.usenonlcon = true;
end

pSc00 = pSc0( : );

%% Determine if need to run opt
constraintTolerance = 1e-9;
constfeasible = checkConstFeasible(...
    pSc00, opt, efsSettings.usenonlcon, constraintTolerance );

%% Loop the ensure feasible start optimization
if ~all( constfeasible )
    startTime = datetime;

    fprintf( "\n----------------------------------------------------------------\n" )
    fprintf( "Starting ensureFeasibleStartFMIN:\t%s", string(datetime) )
    fprintf( "\n----------------------------------------------------------------\n" )

    ii = 1;
    maxIter = 5;

    returnSamePtTrcker = 0;
    returnSamePtTrckerMax = 3;

    maxLoopOptTime = 10 * 60;

    maxEndScale = 0.005;

    while true

        timeLeft = efsSettings.maxTime - seconds( datetime - startTime );

        fprintf( "\n-------------------------\n" );
        fprintf( "iter:\t%i\n", ii );
        fprintf( "time:\t%s\n", string(datetime('now')) );
        fprintf( "time left: %s\n", string( round( timeLeft ) ) );

        timeConst = min( [ timeLeft, maxLoopOptTime ] );

        % [ pSc0Opt, success, returnSamePt, feasStartOutput ] = loopOpt(...
        %     pSc0, pSc00, opt, efsSettings.usenonlcon, timeConst, efsSettings.relLineSrchBnd );

        [ pSc0Opt, success, returnSamePt, feasStartOutput ] = loopOpt(...
            pSc0, pSc0, opt, efsSettings.usenonlcon, timeConst, efsSettings.relLineSrchBnd );

        fprintf( "\ntotal x0 deviation: %s\n", string( norm(pSc0Opt - pSc00 )/norm(pSc00) ) );
        fprintf( "\n-------------------------\n" );

        currElapsedTime = seconds( datetime - startTime );

        if returnSamePt
            efsSettings.relLineSrchBnd = 2.0 * efsSettings.relLineSrchBnd;
            returnSamePtTrcker = returnSamePtTrcker + 1;
        end

        if currElapsedTime > efsSettings.maxTime
            warning( "ensureFeasibleStartFMIN.m could not find feasible point." )
            timeExceed = true;
        else
            timeExceed = false;
        end

        if ( success )
            break;
        end

        if timeExceed || ( ii >= maxIter ) || ( returnSamePtTrcker >= returnSamePtTrckerMax )
            
            if efsSettings.scaleAtEnd
                [ ~, maxConstViolEnd ] = checkConstFeasible(...
                    pSc0Opt, opt, efsSettings.usenonlcon, constraintTolerance );
                pSc0Opt = ( 1 - min( [ maxConstViolEnd, maxEndScale ] ) ) * pSc0Opt;
            end
            break;
        end

        pSc0 = pSc0Opt;
        ii = ii + 1;
    end

    fprintf( "\n----------------------------------------------------------------\n" )
    fprintf( "Finishing ensureFeasibleStartFMIN:\t%s", string(datetime) )
    fprintf( "\n----------------------------------------------------------------\n" )
else
    pSc0Opt = pSc00;
    feasStartOutput.fval = 0;
end


end

%% Helper Functions
% ----------------------------------------------------------------------- %
function [ pSc0Opt, stop, returnSamePt, feasStartOutput ] = loopOpt(...
    pSc0, pSc00, opt, usenonlcon, timeConst, relLineSrchBnd )

%% Prepare Opt
pSc0 = pSc0( : );
clear costFnTrackDecisionVariablesWrapper;
clear feasStarFun;
feasStarFun = @( pSci ) feasStartObjFun( pSci, double( pSc00 ) );

%% Run Opt
A = full( opt.A );
b = opt.b;
Aeq = full( opt.Aeq );
beq = opt.beq;
lb = opt.lb;
ub = opt.ub;
if usenonlcon
    nonlcon = opt.nonlcon;
else
    nonlcon = [];
end

fminopt = optimoptions( "fmincon" );
fminopt.SpecifyConstraintGradient = true;
fminopt.SpecifyObjectiveGradient = true;
fminopt.MaxFunctionEvaluations = inf;
fminopt.FunctionTolerance = 1e-5;
fminopt.OptimalityTolerance = 1e-7;
fminopt.Display = "iter";

fminopt.ConstraintTolerance = 1e-9;

fminopt.MaxIterations = 1000;

% fminopt.Algorithm = "sqp";

fminopt.StepTolerance = 1e-6;
fminopt.Algorithm = "active-set";
fminopt.RelLineSrchBnd = relLineSrchBnd;
fminopt.RelLineSrchBndDuration = inf;
fminopt.TolConSQP = 1e-12;

% fminopt.Algorithm = "interior-point";
% fminopt.StepTolerance = 1e-10;
% fminopt.EnableFeasibilityMode = true;
% fminopt.SubproblemAlgorithm = 'cg';
% fminopt.BarrierParamUpdate = 'predictor-corrector';
% fminopt.InitBarrierParam = 1e-2;
% fminopt.TolProjCGAbs = 1e-12;
% fminopt.TolProjCG = 1e-8;

fminopt.OutputFcn = prepareMaxTimeFMin( timeConst );

feasStartOutput = struct;

pSc0 = double( pSc0 );

[ constfeasible, maxConstViolSt ] = checkConstFeasible(...
    pSc0, opt, usenonlcon, fminopt.ConstraintTolerance );

if all( constfeasible )
    clear costFnTrackDecisionVariablesWrapper;
    clear feasStarFun;
    
    stop = true;
    feasStartOutput.fval = 0;
    pSc0Opt = pSc0;
else

    spacingTrackDecisionVariables = 1;
    constrTolSave = 1e3;

    const = struct;
    const.spacingTrackDecisionVariables = spacingTrackDecisionVariables;
    const.constrTolSave = constrTolSave;
    const.A = A;
    const.b = b;
    const.Aeq = Aeq;
    const.beq = beq;
    const.nonlcon = nonlcon;
    
    feasStarFun = @( pSc ) costFnTrackDecisionVariablesWrapper(...
        pSc, feasStarFun, const, const.spacingTrackDecisionVariables, const.constrTolSave );

    [ pSc0Opt, feasStartOutput.fval, feasStartOutput.exitflag, feasStartOutput.output ] =...
        fmincon( feasStarFun, pSc0,...
        A, b, Aeq, beq, lb, ub, nonlcon, fminopt );

    [ ~, ~, trackDecSt ] =...
        feasStarFun( pSc0Opt );
    clear costFnTrackDecisionVariablesWrapper;
    clear feasStarFun;

    if ~isempty( feasStartOutput.output.bestfeasible )
        pSc0Opt = feasStartOutput.output.bestfeasible.x( : );
        feasStartOutput.fval = feasStartOutput.output.bestfeasible.fval;
        stop = true;
        returnSamePt = false;
    else
        [ ~, minConstViolIdx ] = min( trackDecSt.const_vars_iters );
        pSc0Opt = trackDecSt.xtrack_vars_iters( :, minConstViolIdx );
        
        [ ~, maxConstViolEnd ] = checkConstFeasible( pSc0Opt, opt, usenonlcon, fminopt.ConstraintTolerance );

        stop = false;

        relConstCh = abs( ( maxConstViolEnd - maxConstViolSt ) / maxConstViolSt );

        if relConstCh < 1e-8
            returnSamePt = true;
        else
            returnSamePt = false;
        end

    end
    

end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ constfeasible, maxConstViol ] = checkConstFeasible(...
    pSc0, opt, usenonlcon, constraintTolerance )

if nargout > 1
    maxConstViol = 0;
end

A = full( opt.A );
b = opt.b;
Aeq = full( opt.Aeq );
beq = opt.beq;
lb = opt.lb;
ub = opt.ub;
if usenonlcon
    nonlcon = opt.nonlcon;
else
    nonlcon = [];
end

pSc0 = double( pSc0 );

% Check to see if constraints are feasible before running the optimization
if ~isempty( A )
    AbEval = ( A * pSc0 - b );
    Abfeasible = all( AbEval <= constraintTolerance );
    if nargout > 1
        if ~Abfeasible
            maxConstViol = max( [ maxConstViol; AbEval ] );
        end
    end
else
    Abfeasible = true;
end
if ~isempty( Aeq )
    AbeqEval = abs( Aeq * pSc0 - beq );
    Abeqfeasible = all( abs( Aeq * pSc0 - beq ) <= constraintTolerance );
    if nargout > 1
        if ~Abeqfeasible
            maxConstViol = max( [ maxConstViol; AbeqEval ] );
        end
    end
else
    Abeqfeasible = true;
end

lbEval = ( lb - pSc0 );
lbfeasible = all( lbEval <= constraintTolerance );

if nargout > 1
    if ~lbfeasible
        maxConstViol = max( [ maxConstViol; lbEval ] );
    end
end

ubEval = ( pSc0 - ub );
ubfeasible = all( ubEval <= constraintTolerance );

if nargout > 1
    if ~ubfeasible
        maxConstViol = max( [ maxConstViol; ubEval ] );
    end
end

if ~isempty( nonlcon )
    [ c, ceq ] = nonlcon( pSc0 );

    if ~isempty( c )
        nonlconIneqfeasible = all( c <= constraintTolerance );
        
        if nargout > 1
            if ~nonlconIneqfeasible
                maxConstViol = max( [ maxConstViol; c ] );
            end
        end

    else
        nonlconIneqfeasible = true;
    end
    if ~isempty( ceq )
        absceq = abs( ceq );
        nonlconEqfeasible = all( absceq <= constraintTolerance );
        
        if nargout > 1
            if ~nonlconEqfeasible
                maxConstViol = max( [ maxConstViol; absceq ] );
            end
        end

    else
        nonlconEqfeasible = true;
    end

else
    nonlconIneqfeasible = true;
    nonlconEqfeasible = true;
end

constfeasible = [...
    Abfeasible;...
    Abeqfeasible;...
    lbfeasible;...
    ubfeasible;...
    nonlconIneqfeasible;...
    nonlconEqfeasible;...
    ];

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ cost, gradcost ] = feasStartObjFun( pSci, pSc0 )
dpSc = pSci - pSc0;

denom = norm( pSc0 );

if denom ~= 0
    cost = norm( dpSc ) / denom;
else
    cost = norm( dpSc );
end

if nargout > 1
    if any( dpSc ~= 0 )
        gradcost = dpSc / ( norm(dpSc) * denom );
    elseif denom == 0
        gradcost = dpSc / ( norm(dpSc) );
    elseif all( dpSc == 0 )
        gradcost = ones( size( pSci ) ) / denom;
    else
        error( "Don't know how to process cost gradient." )
    end
end
end
% ----------------------------------------------------------------------- %



% Old code
%
%
% function [ pSc0Opt, feasStartOutput ] = ensureFeasibleStartFMIN( pSc0, opt, usenonlcon, maxTime, relLineSrchBnd )
% if nargin < 3
%     usenonlcon = true;
% end
% if nargin < 4
%     maxTime = 30 * 60; % seconds
% end
%
% if nargin < 5
%   % relLineSrchBnd = max( min( 1e-4 * sqrt( double( opt.numVars ) / 250 ), 1e-3 ), 5e-4 );
%   % relLineSrchBnd = 1e-4;
%   % relLineSrchBnd = 8.75e-5;
%   relLineSrchBnd = 5e-5;
% end
%
% %% Prepare Opt
% pSc0 = pSc0( : );
% feasStarFun = @( pSci ) feasStartObjFun( pSci, double( pSc0 ) );
%
% %% Run Opt
% A = full( opt.A );
% b = opt.b;
% Aeq = full( opt.Aeq );
% beq = opt.beq;
% lb = opt.lb;
% ub = opt.ub;
% if usenonlcon
%     nonlcon = opt.nonlcon;
% else
%     nonlcon = [];
% end
%
% fminopt = optimoptions( "fmincon" );
% fminopt.SpecifyConstraintGradient = true;
% fminopt.SpecifyObjectiveGradient = true;
% fminopt.MaxFunctionEvaluations = inf;
% fminopt.ConstraintTolerance = 1e-9;
% fminopt.FunctionTolerance = 1e-5;
% fminopt.OptimalityTolerance = 1e-7;
% fminopt.Display = "iter";
%
% fminopt.MaxIterations = 2e3;
%
% % fminopt.Algorithm = "sqp";
%
% fminopt.StepTolerance = 1e-6;
% fminopt.Algorithm = "active-set";
% fminopt.RelLineSrchBnd = relLineSrchBnd;
% fminopt.RelLineSrchBndDuration = inf;
% fminopt.TolConSQP = 1e-8;
% % fminopt.MaxSQPIter = 5e3;
%
% fminopt.OutputFcn = prepareMaxTimeFMin( maxTime );
%
% % if matches( opt.fminopt.Algorithm, "interior-point", 'IgnoreCase', true )
% %
% %     fminopt.Algorithm = "interior-point";
% %     fminopt.StepTolerance = 1e-10;
% %     fminopt.EnableFeasibilityMode = true;
% %     fminopt.SubproblemAlgorithm = 'cg';
% %     fminopt.BarrierParamUpdate = 'predictor-corrector';
% %     fminopt.InitBarrierParam = 1e-1;
% %     fminopt.TolProjCGAbs = 1e-12;
% %     fminopt.TolProjCG = 1e-4;
% %
% % elseif matches( opt.fminopt.Algorithm, [ "active-set"; "sqp" ], 'IgnoreCase', true )
% %
% %     % fminopt.Algorithm = "sqp";
% %
% %     fminopt.StepTolerance = 1e-6;
% %     fminopt.Algorithm = "active-set";
% %     fminopt.RelLineSrchBnd = 1e-4;
% %     fminopt.RelLineSrchBndDuration = inf;
% %     fminopt.TolConSQP = 1e-12;
% % end
%
% feasStartOutput = struct;
%
% pSc0 = double( pSc0 );
%
% % Check to see if constraints are feasible before running the optimization
% if ~isempty( A )
%     Abfeasible = all( ( A * pSc0 - b ) <= fminopt.ConstraintTolerance );
% else
%     Abfeasible = true;
% end
% if ~isempty( Aeq )
%     Abeqfeasible = all( abs( Aeq * pSc0 - beq ) <= fminopt.ConstraintTolerance );
% else
%     Abeqfeasible = true;
% end
% lbfeasible = all( ( pSc0 - lb ) >= fminopt.ConstraintTolerance );
% ubfeasible = all( ( pSc0 - ub ) <= fminopt.ConstraintTolerance );
% if ~isempty( nonlcon )
%     [ c, ceq ] = nonlcon( pSc0 );
%     if ~isempty( c )
%         nonlconIneqfeasible = all( c <= fminopt.ConstraintTolerance );
%     else
%         nonlconIneqfeasible = true;
%     end
%     if ~isempty( ceq )
%         nonlconEqfeasible = all( abs( ceq ) <= fminopt.ConstraintTolerance );
%     else
%         nonlconEqfeasible = true;
%     end
%
% else
%     nonlconIneqfeasible = true;
%     nonlconEqfeasible = true;
% end
% constfeasible = [ Abfeasible; Abeqfeasible; lbfeasible; ubfeasible; nonlconIneqfeasible; nonlconEqfeasible ];
%
% if all( constfeasible )
%     clear costFnTrackDecisionVariablesWrapper;
%     clear feasStarFun;
%
%     feasStartOutput.fval = 0;
%     pSc0Opt = pSc0;
% else
%
%     fprintf( "\n----------------------------------------------------------------\n" )
%     fprintf( "Starting ensureFeasibleStartFMIN:\t%s", string(datetime) )
%     fprintf( "\n----------------------------------------------------------------\n" )
%
%     spacingTrackDecisionVariables = 1;
%     constrTolSave = 5e-2;
%
%     const = struct;
%     const.spacingTrackDecisionVariables = spacingTrackDecisionVariables;
%     const.constrTolSave = constrTolSave;
%     const.A = A;
%     const.b = b;
%     const.Aeq = Aeq;
%     const.beq = beq;
%     const.nonlcon = nonlcon;
%
%     feasStarFun = @( pSc ) costFnTrackDecisionVariablesWrapper(...
%         pSc, feasStarFun, const, const.spacingTrackDecisionVariables, const.constrTolSave );
%
%     [ pSc0Opt, feasStartOutput.fval, feasStartOutput.exitflag, feasStartOutput.output ] =...
%         fmincon( feasStarFun, pSc0,...
%         A, b, Aeq, beq, lb, ub, nonlcon, fminopt );
%
%     [ ~, ~, trackDecSt ] =...
%         feasStarFun( pSc0Opt );
%     clear costFnTrackDecisionVariablesWrapper;
%     clear feasStarFun;
%
%     fprintf( "\n----------------------------------------------------------------\n" )
%     fprintf( "Finishing ensureFeasibleStartFMIN:\t%s", string(datetime) )
%     fprintf( "\n----------------------------------------------------------------\n" )
%
%     if ~isempty( feasStartOutput.output.bestfeasible )
%         pSc0Opt = feasStartOutput.output.bestfeasible.x( : );
%         feasStartOutput.fval = feasStartOutput.output.bestfeasible.fval;
%     else
%         [ ~, minConstViolIdx ] = min( trackDecSt.const_vars_iters );
%         pSc0Opt = trackDecSt.xtrack_vars_iters( :, minConstViolIdx );
%         warning( "ensureFeasibleStartFMIN.m could not find feasible point." )
%     end
% end
%
% end
%
% %% Helper Functions
% % ----------------------------------------------------------------------- %
% function [ cost, gradcost ] = feasStartObjFun( pSci, pSc0 )
% dpSc = pSci - pSc0;
%
% denom = norm( pSc0 );
%
% if denom ~= 0
%     cost = norm( dpSc ) / denom;
% else
%     cost = norm( dpSc );
% end
%
% if nargout > 1
%     if any( dpSc ~= 0 )
%         gradcost = dpSc / ( norm(dpSc) * denom );
%     elseif denom == 0
%         gradcost = dpSc / ( norm(dpSc) );
%     elseif all( dpSc == 0 )
%         gradcost = ones( size( pSci ) ) / denom;
%     else
%         error( "Don't know how to process cost gradient." )
%     end
% end
% end
% % ----------------------------------------------------------------------- %