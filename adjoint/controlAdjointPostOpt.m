function [ opt, valtrain, valtest, oc, pulse ] = controlAdjointPostOpt(...
    opt, oc, pulse, fields )
% This function will control the adjoint optimization process. It will
% first process the pulses and constraints. It will then determine the
% number of optimization variables. It will then generate the train and
% test databases that might be used for universal pulse design.

%% Process fields for val struct
timeResDigits = 6;
valtrain = struct;
valtrain.structtype = 'val';
valtrain.di = oc.val_di;
valtrain.dt = round( oc.val_dt, timeResDigits);
% valtrain.opt_tvec = oc.opt_tvec;

[ oc, ~, ~, valtrain ] = processFieldsStruct(...
    oc, pulse, fields, fields.si_train, valtrain );
[ opt, valtrain ] = addOptInfoVal( opt, valtrain );
[ oc, pulse, valtrain ] = processPulseName( oc, pulse, valtrain );
[ oc, ~, valtrain ] = processAdjointFunctions( oc, pulse, valtrain );
valtrain.si = fields.si_train;
valtrain.flipAngleRangePlot = pulse.flipAngleRangePlot;

if ~isempty( fields.si_test )
    valtest = struct;
    valtest.structtype = 'val';
    valtest.di = oc.val_di;
    valtest.dt = valtrain.dt;
    % valtest.opt_tvec = oc.opt_tvec;
    valtest.si = fields.si_test;
    valtest.flipAngleRangePlot = pulse.flipAngleRangePlot;
    
    [ oc, ~, ~, valtest ] = processFieldsStruct(...
        oc, pulse, fields, fields.si_test, valtest );
    [ opt, valtest ] = addOptInfoVal( opt, valtest );
    [ oc, pulse, valtest ] = processPulseName( oc, pulse, valtest );
    [ oc, ~, valtest ] = processAdjointFunctions( oc, pulse, valtest );

else
    valtest = [];
end

%% Post process
[ valtrain, opt ] = opt.postProcessFunction( valtrain, opt, pulse, fields );

if ~isempty( valtest )
    [ valtest, opt ] = opt.postProcessFunction( valtest, opt, pulse, fields );
end

clear fields;

%% Process Adjoint Images
[ cellFigs, strLabels ] = processAdjointImages( opt, pulse, valtrain, valtest, oc ); %#ok

%% Process Adjoint Data
if oc.saveResult
    note = collateAdjointData( valtrain, valtest, opt, pulse, oc );
    if isfield( oc, 'currFile' ) && isfield( oc, 'saveFileRecord' ) && oc.saveFileRecord
        fileRecord = getFunctionRecord( oc.currFile );
    else
        oc.saveFileRecord = false;
    end
end

%% Clean up opt struct
if oc.saveResult
    
    optPreserve = opt;
    valtrainPreserve = valtrain;
    pulsePreserve = pulse;

    opt = cleanUpAdjointStructs( optPreserve );
    valtrain = cleanUpAdjointStructs( valtrainPreserve );

    if ~isempty( valtest )
        valtestPreserve = valtest;
        valtest = cleanUpAdjointStructs( valtestPreserve );
    end

    % Remove M0 and Mtarg for easy saving
    if isfield( pulse, 'M0' )
        pulse = rmfield( pulse, 'M0' );
    end
    if isfield( pulse, 'Mtarg' )
        pulse = rmfield( pulse, 'Mtarg' );
    end
    if isfield( pulse, 'targSt' )
        if isfield( pulse.targSt, 'MtargVals' )
            if numel( pulse.targSt.MtargVals ) > 3
                pulse.targSt = rmfield( pulse.targSt, 'MtargVals' );
            end
        end
    end

end

%% Save Optimization Results
if oc.saveResult

    % save result
    % Change warning temporarily
    toobigWarn = warning('error', 'MATLAB:save:sizeTooBigForMATFile');

    if ~isempty( valtest )
        try
            save( fullfile( oc.saveDir, 'optOutput.mat' ),...
                'opt', 'pulse', 'oc', 'valtrain', 'valtest' );
        catch
            save( fullfile( oc.saveDir, 'optOutput.mat' ),...
                'opt', 'pulse', 'oc', 'valtrain', 'valtest', "-v7.3" );
        end
    else
        try
            save( fullfile( oc.saveDir, 'optOutput.mat' ),...
                'opt', 'pulse', 'oc', 'valtrain' );
        catch
            save( fullfile( oc.saveDir, 'optOutput.mat' ),...
                'opt', 'pulse', 'oc', 'valtrain', "-v7.3" );
        end
    end

    % Change warning back
    toobigWarn.state = 'on';

    % save notes
    saveTxtFromStrArray( fullfile( oc.saveDir, 'note.txt' ), note );

    if oc.saveFileRecord
        fileRecord = [ strcat( "% ", oc.currFile ); ""; ""; fileRecord ];
        saveTxtFromStrArray( fullfile( oc.saveDir, 'fileRecord.txt' ), fileRecord );
    end
end

%% Regenerate structs
if oc.saveResult
    opt = optPreserve;
    clear optPreserve;

    valtrain = valtrainPreserve;
    clear valtrainPreserve;

    if ~isempty( valtest )
        valtest = valtestPreserve;
        clear valtestPreserve;
    end

    pulse = pulsePreserve;
    clear pulsePreserve;
end

end