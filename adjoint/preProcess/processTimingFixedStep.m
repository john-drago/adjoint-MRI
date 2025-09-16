function [ opt, oc, tvec, dtvec, numTimePoints, pulseLength ] =...
    processTimingFixedStep( oc, pulse, opt )


timeResDigits = 6; % round to microseconds
pulseLength = round( pulse.length, timeResDigits); 
dt = opt.dt; % get spacing of points for optimization
dt_tol = dt / 10;

% Determine number of full periods
numFullPeriods = floor( ( pulseLength - dt_tol ) / dt );

tvec = (dt/2) : dt : ( numFullPeriods * dt + dt_tol );
dtvec = dt * ones( size( tvec ) );

dtLast = ( pulseLength - numFullPeriods * dt );

tvec = [ tvec, numFullPeriods * dt + dtLast/2 ];
dtvec = [ dtvec, dtLast ];

tvec = round( tvec, timeResDigits+1 );
dtvec = round( dtvec, timeResDigits+1 );

numTimePoints = length( tvec );

opt.numTimePoints = numTimePoints;
opt.tvec = tvec;
opt.dtvec = dtvec;
opt.pulseLength = pulseLength;

if strcmpi( opt.structtype, "val" )
    opt.opt_tvec = oc.opt_tvec;
    opt.opt_dtvec = oc.opt_dtvec;
elseif strcmpi( opt.structtype, "opt" )
    opt.opt_tvec = tvec;
    opt.opt_dtvec = dtvec;
    oc.opt_tvec = tvec;
    oc.opt_dtvec = dtvec;
end

end