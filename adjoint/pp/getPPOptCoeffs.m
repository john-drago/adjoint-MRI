function c = getPPOptCoeffs( wv, opt, tol )

if nargin < 3
    tol = 1e-10;
end

if ~isfield( wv, 'breal' ) || ~isfield( wv, 'bimag' )
    wv.breal = real( wv.RF );
    wv.bimag = imag( wv.RF );
end

%% Process RF waveforms
% Prepare struct to pass in
RFSt = getProjSt( opt, 'RF' );

c.breal_coeffs = getPPCoeffWaveform( wv.tvec, wv.breal, RFSt, tol );
c.bimag_coeffs = getPPCoeffWaveform( wv.tvec, wv.bimag, RFSt, tol );

%% Process Grad waveforms
% Prepare struct to pass in
gradSt = getProjSt( opt, 'grad' );

c.grad_coeffs = getPPCoeffWaveform( wv.tvec, wv.Grad, gradSt, tol );

%% Process Shim waveforms
if opt.numZCoils > 0
    % Prepare struct to pass in
    shimSt = getProjSt( opt, 'shim' );

    c.shim_coeffs = getPPCoeffWaveform( wv.tvec, wv.Shim, shimSt, tol );
else
    c.shim_coeffs = zeros( 1, opt.numZCoils );
end

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function projSt = getProjSt( opt, label )

projSt = struct;
projSt.numPPPeriods = opt.( sprintf( "numPPPeriods_%s", label ) );
projSt.PPStartStop = opt.( sprintf( "PPStartStop_%s", label ) );
projSt.PPVarIdxs = opt.( sprintf( "PPVarIdxs_%s", label ) );
projSt.numVarsPerChannel = opt.( sprintf( "numVarsPerChannel_%s", label ) );
projSt.numPPShape = opt.( sprintf( "numPPShape_%s", label ) );
projSt.PPVarSFIdxs = opt.( sprintf( "PPVarSFIdxs_%s", label ) );
projSt.PPSFIdxs = opt.( sprintf( "PPSFIdxs_%s", label ) );
projSt.PPVarIdxsCell = opt.( sprintf( "PPVarIdxsCell_%s", label ) );
projSt.PPVarIdxsCoilCell = opt.( sprintf( "PPVarIdxsCoilCell_%s", label ) );
projSt.shapeFnChebCoeffs = opt.( sprintf( "shapeFnChebCoeffs_%s", label ) );
projSt.PPIdxs = opt.( sprintf( "PPIdxs_%s", label ) );
projSt.orderPP = opt.( sprintf( "orderPP_%s", label ) );
projSt.varsToTimepoints = opt.( sprintf( "varsToTimepoints_%s", label ) );
projSt.tvec = opt.tvec;

end
% ----------------------------------------------------------------------- %