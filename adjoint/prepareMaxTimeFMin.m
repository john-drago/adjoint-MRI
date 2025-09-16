function OutputFcn = prepareMaxTimeFMin( maxTime )

OutputFcn = @( x, optimValues, state )...
    checkMaxTimeFMin( x, optimValues, state, maxTime );

end