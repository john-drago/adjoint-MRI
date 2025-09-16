function [ oc, opt, pulse, pulseopt, fields, valtrain, valtest ] = ...
    recreateOptStructs( optSave )
% Function will recreate the structs that are used during optimization from
% the output of the optOutput.mat file.

%% Get the structs
oc = optSave.oc;
pulse = optSave.pulse;
opt = optSave.opt;
valtrain = optSave.valtrain;
if isfield( optSave, 'valtest' )
    valtest = optSave.valtest;
else
    valtest = [];
end

%% Start wih field data
fields = struct;

% subject maps
subjectMapsPath = regenerateFilePathOpt( pulse.subjectMaps );
fields.subjectMaps = subjectMapsPath;

% process UP
UP = load( subjectMapsPath );

% determine if tailored or not
if contains( pulse.optPopulation, "tailored", 'IgnoreCase', true )
    if isfield( pulse, "tailoredTrain" )
        idenLogical = contains( UP.subjIden, pulse.tailoredTrain, 'ignorecase', true );
    else
        idenLogical = true;
    end
elseif contains( pulse.optPopulation, "universal", 'IgnoreCase', true )
    idenLogical = contains( UP.subjIden, [ pulse.universalTrain; pulse.universalTest ],...
        'ignorecase', true );
end

fields.subjIden = pulse.subjIden;
fields.subjPath = pulse.subjPath;

fields.x = UP.x;
fields.y = UP.y;
fields.z = UP.z;
fields.X = UP.X;
fields.Y = UP.Y;
fields.Z = UP.Z;
if opt.numXYCoils > 1
    fields.b1p = UP.b1p( :, :, :, :, idenLogical );
else
    fields.b1p = UP.b1p( :, :, :, idenLogical );
end
fields.db0 = UP.db0shim( :, :, :, idenLogical );
fields.opt_roi = UP.roi_brain( :, :, :, idenLogical );
fields.val_roi = UP.roi_body( :, :, :, idenLogical );

% VOPs
if isfield( pulse, 'VOPpath' ) && ~isempty( pulse.VOPpath )
    VOPpath = regenerateFilePathOpt( pulse.VOPpath );
    fields.VOPpath = VOPpath;

    VOPst = load( VOPpath );
    fields.VOPs = VOPst.VOPs;

    if isfield( VOPst, 'QPartialBody' )
        fields.QGlobal = VOPst.QPartialBody;
    elseif isfield( VOPst, 'QGlobalMat' )
        fields.QGlobal = VOPst.QGlobalMat;
    end
end

% shim coil maps
if isfield( pulse, 'zCoilPath' ) && ~isempty( pulse.zCoilPath )
    zCoilPath = regenerateFilePathOpt( pulse.zCoilPath );
    fields.zCoilPath = zCoilPath;

    BS = load( zCoilPath );

    if isfield( pulse, 'zCoilShiftPath' ) && ~isempty( pulse.zCoilShiftPath )
        zCoilShiftPath = regenerateFilePathOpt( pulse.zCoilShiftPath );
        fields.zCoilShiftPath = zCoilShiftPath;
        shift = load( zCoilShiftPath );

        BS.X = BS.X + shift.shiftArray( 1, 1 );
        BS.Y = BS.Y + shift.shiftArray( 1, 2 );
        BS.Z = BS.Z + shift.shiftArray( 1, 3 );
        BS = rmfield( BS, ["x", "y", "z"] );
    end


    fields.bz = zeros( [ BS.coilNum, size( fields.X ) ] );

    for cc = 1:BS.coilNum
        fields.bz( cc, :, :, : ) = Interp3D(...
            BS.X, BS.Y, BS.Z, squeeze(BS.bz( cc, :, :, : )),...
            UP.X, UP.Y, UP.Z );
    end

else
    fields.bz = 0;
end

%% Process fields struct 
% Get data about the fields such as number of coils and number of subjects
[ oc, pulse, fields ] = processFieldsData( oc, pulse, fields );

% Now add fields to the opt struct
[ oc, pulse, fields, opt ] = processFieldsStruct(...
    oc, pulse, fields, fields.si_train, opt );

%% Deal with opt struct
[ oc, pulseopt, opt ] = processPulseName( oc, pulse, opt );

%% Add functions to structs
[ oc, pulseopt, opt ] = processAdjointFunctions( oc, pulseopt, opt );

ensureFeasibleStart = false;
[ oc, pulseopt, opt ] = prepareAdjointOpt(...
    oc, pulseopt, opt, ensureFeasibleStart );

%% Deal with valtrain
[ oc, ~, ~, valtrain ] = processFieldsStruct(...
    oc, pulseopt, fields, fields.si_train, valtrain );
[ opt, valtrain ] = addOptInfoVal( opt, valtrain );
[ oc, pulsevaltrain, valtrain ] = processPulseName( oc, pulse, valtrain );
[ oc, ~, valtrain ] = processAdjointFunctions( oc, pulsevaltrain, valtrain );
valtrain.si_train = fields.si_train;
valtrain.flipAngleRangePlot = pulseopt.flipAngleRangePlot;

[ valtrain, opt ] = opt.postProcessFunction( valtrain, opt, pulse );

%% Deal with valtest
if ~isempty( valtest )
    [ oc, ~, ~, valtest ] = processFieldsStruct(...
        oc, pulseopt, fields, fields.si_test, valtest );
    [ opt, valtest ] = addOptInfoVal( opt, valtest );
    [ oc, pulsevaltest, valtest ] = processPulseName( oc, pulse, valtest );
    [ oc, ~, valtest ] = processAdjointFunctions( oc, pulsevaltest, valtest );
    valtest.si_test = fields.si_test;
    valtest.flipAngleRangePlot = pulseopt.flipAngleRangePlot;

    [ valtest, opt ] = opt.postProcessFunction( valtest, opt, pulse );
else
    valtest = [];
end

end