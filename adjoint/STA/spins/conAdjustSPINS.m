function [ AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt ] = conAdjustSPINS(...
    opt, staSt, AbSt, AbeqSt, nlconIneqSt, nlconEqSt, binInit )
% This file will add in nonlinear constraints for the SPINS pulses to
% ensure compliance

%% Determine if there is a grad max
if isfield( opt, 'gradMax_constr' )
    [ AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt ] = ...
            determineSPINSConstraintFunction(...
            'grad-max', AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt, opt, binInit ); 
end

%% Iterate over nonlinear inequality constraints first
if ~isempty( opt.nlconIneqFuncs )
    for nn = 1:length( opt.nlconIneqFuncs )

        spinsConString = func2str( opt.nlconIneqFuncs{ nn, 1 } );
        [ AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt ] = ...
            determineSPINSConstraintFunction(...
            spinsConString, AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt, opt, binInit ); 

    end
end

%% Iterate over nonlinear equality constraints
if ~isempty( opt.nlconEqFuncs )
    for nn = 1:length( opt.nlconEqFuncs )

        spinsConString = func2str( opt.nlconEqFuncs{ nn, 1 } );
        [ AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt ] = ...
            determineSPINSConstraintFunction(...
            spinsConString, AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt, opt, binInit ); 

    end
end

%% Iterate over linear inequality constraints
if ~isempty( opt.AbConstraintNames )
    for nn = 1:length( opt.AbConstraintNames )
        
        [ AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt ] = ...
            determineSPINSConstraintFunction(...
            opt.AbConstraintNames( nn ), AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt, opt, binInit ); 

    end
end

%% Iterate over linear equality constraints
if ~isempty( opt.AbeqConstraintNames )
    for nn = 1:length( opt.AbeqConstraintNames )

        [ AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt ] = ...
            determineSPINSConstraintFunction(...
            opt.AbeqConstraintNames( nn ), AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt, opt, binInit );  

    end
end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function [ AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt ] = determineSPINSConstraintFunction(...
    spinsConString, AbSt, AbeqSt, nlconIneqSt, nlconEqSt, staSt, opt, binInit )

if contains( spinsConString, [ "TotalRFPower"; "total-RF-power" ], 'ignorecase', true )
    
    spinsConstraintFunc = @constraintSPINSPulseTotalRFPower;
    staSt.totalRFPower_constr = opt.totalRFPower_constr;

    [ nlconIneqSt ] = nlconAddList(...
                nlconIneqSt, spinsConstraintFunc, 1 );

elseif contains( spinsConString, [ "MaxRFPower"; "max-RF-power" ], 'ignorecase', true )

    spinsConstraintFunc = @constraintSPINSPulseMaxRFPower;
    staSt.maxRFPower_constr = opt.maxRFPower_constr;

    [ nlconIneqSt ] = nlconAddList(...
                nlconIneqSt, spinsConstraintFunc, 1 );

elseif contains( spinsConString, [ "PeakLocalSAR"; "peak-local-SAR" ], 'ignorecase', true )

    spinsConstraintFunc = @constraintSPINSPulsePeakLocalSAR;
    staSt.peakLocalSAR_constr = opt.peakLocalSAR_constr;

    [ nlconIneqSt ] = nlconAddList(...
                nlconIneqSt, spinsConstraintFunc, 1 );

elseif contains( spinsConString, [ "PeakGlobalSAR"; "peak-global-SAR" ], 'ignorecase', true )

    spinsConstraintFunc = @constraintSPINSPulsePeakGlobalSAR;
    staSt.peakGlobalSAR_constr = opt.peakGlobalSAR_constr;

    [ nlconIneqSt ] = nlconAddList(...
                nlconIneqSt, spinsConstraintFunc, 1 );

elseif contains( spinsConString, [ "AvgLocalSAR"; "average-local-SAR" ], 'ignorecase', true )

    spinsConstraintFunc = @constraintSPINSPulseAvgLocalSAR;
    staSt.avgLocalSAR_constr = opt.avgLocalSAR_constr;

    [ nlconIneqSt ] = nlconAddList(...
                nlconIneqSt, spinsConstraintFunc, opt.numVOPs );

elseif contains( spinsConString, [ "AvgGlobalSAR"; "average-global-SAR" ], 'ignorecase', true )

    spinsConstraintFunc = @constraintSPINSPulseAvgGlobalSAR;
    staSt.avgGlobalSAR_constr = opt.avgGlobalSAR_constr;

    [ nlconIneqSt ] = nlconAddList(...
                nlconIneqSt, spinsConstraintFunc, 1 );

elseif contains( spinsConString, [ "GradMax"; "grad-max" ], 'ignorecase', true )
    
    if ~binInit
        spinsConstraintFunc = @constraintSPINSPulseGradMax;
        staSt.gradMax_constr = opt.gradMax_constr;

        [ nlconIneqSt ] = nlconAddList(...
                nlconIneqSt, spinsConstraintFunc, 3 );
    end

elseif contains( spinsConString, [ "GradSlewRate"; "grad-slew-rate" ], 'ignorecase', true )
    
    if ~binInit
        spinsConstraintFunc = @constraintSPINSPulseGradSlewRate;
        staSt.gradSlewRate_constr = opt.gradSlewRate_constr;

        [ nlconIneqSt ] = nlconAddList(...
                nlconIneqSt, spinsConstraintFunc, 3 * staSt.spins_int_num + 2 * 3 );
    end

elseif contains( spinsConString, "RF-slew-rate", 'ignorecase', true )
    
    staSt.RFSlewRate_constr = opt.RFSlewRate_constr;
    [ ARFSlew, bRFSlew ] = constraintSPINSPulseRFSlew(...
                staSt, opt.RFSlewRate_constr );

    AbSt = AbAddList( AbSt, "RF-slew-rate", ARFSlew, bRFSlew );

elseif contains( spinsConString, "RF-accel", 'ignorecase', true )
    
    staSt.RFAccel_constr = opt.RFAccel_constr;
    [ ARFAccel, bRFAccel ] = constraintSPINSPulseRFAccel(...
                opt, opt.RFAccel_constr );

    AbSt = AbAddList( AbSt, "RF-accel", ARFAccel, bRFAccel );

end

end
% ----------------------------------------------------------------------- %
