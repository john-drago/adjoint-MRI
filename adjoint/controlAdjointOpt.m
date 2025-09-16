function [ opt, valtrain, valtest, pulse, oc ] = controlAdjointOpt( oc, pulse, fields )
% This function will control the adjoint optimization process. It will
% first process the pulses and constraints. It will then determine the
% number of optimization variables. It will then generate the train and
% test databases that might be used for universal pulse design.

%% Run all pre optimization processes
[ opt, oc, pulse, fields ] = controlAdjointPreOpt(...
    oc, pulse, fields );

%% Run opt
opt = runAdjointOpt( opt, oc );

%% Post process optimization
[ opt, valtrain, valtest, oc, pulse ] = controlAdjointPostOpt(...
    opt, oc, pulse, fields );

end