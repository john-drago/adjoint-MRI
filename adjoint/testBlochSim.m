function [ errBS, wv, wvBS ] = testBlochSim( p, opt )

%% Generate Waveforms
wv = opt.generateWaveforms( p, opt );

%% Run forward model integration
if opt.useGPU
    [ Marray, ~, wv ] = runAdjointForwardModelGPU( wv, opt );
else
    [ Marray, ~, wv ] = runAdjointForwardModelCPU( wv, opt );
end

if opt.useGPU
    MT = Marray( :, :, end );
else
    MT = Marray{ end };
end

%% Test Bloch Sim
wvBS = opt.generatePlotWaveforms( p, opt );
wvBS.bzsens = opt.bzsens;
wvBS.b1preal = opt.b1preal;
wvBS.b1pimag = opt.b1pimag;
wvBS.db0 = opt.db0;
wvBS.pos = opt.pos;
wvBS.M0 = opt.M0;

MBS = blochSimVec( wvBS );

errBS = norm( MT - MBS ) / norm( MBS );

fprintf('Difference Between Rot Matrix and ODE45 Bloch Sim:\t %.5g\n', errBS);

end