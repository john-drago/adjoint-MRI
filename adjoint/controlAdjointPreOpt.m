function [ opt, oc, pulse, fields ] = controlAdjointPreOpt( oc, pulse, fields )
% This function will control the adjoint optimization process. It will
% first process the pulses and constraints. It will then determine the
% number of optimization variables. It will then generate the train and
% test databases that might be used for universal pulse design.

%% Initialize opt and val structs
opt = struct;
opt.structtype = 'opt';

timeResDigits = 6; % round to microseconds
opt.di = oc.opt_di;
opt.dt = round( oc.opt_dt, timeResDigits);

%% Process fields struct 
% Get data about the fields such as number of coils and number of subjects
[ oc, pulse, fields ] = processFieldsData( oc, pulse, fields );

% Now add fields to the opt struct
[ oc, pulse, fields, opt ] = processFieldsStruct(...
    oc, pulse, fields, fields.si_train, opt );

%% Process pulse struct
[ oc, pulse, opt ] = processPulseName( oc, pulse, opt );

%% Add functions to execute forward model and calculate cost and gradients
[ oc, pulse, opt ] = processAdjointFunctions( oc, pulse, opt );

%% Prepare optimization functions
[ oc, pulse, opt ] = prepareAdjointOpt( oc, pulse, opt );

end