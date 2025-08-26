function [ oc, opt ] = assignGPUInformation( oc, opt )

if isfield( oc, 'useGPU' )
    if isfield( oc, "fminopt" )
        fminUseParallel = oc.fminopt.UseParallel;
    else
        fminUseParallel = false;
    end
    if isfield( oc, "gaopt" )
        gaUseParallel = oc.fminopt.UseParallel;
    else
        gaUseParallel = false;
    end
    if ~matches( oc.optType, "ga-proxy", "ignorecase", true )
        if ( fminUseParallel || gaUseParallel ) && oc.useGPU
            opt.useGPU = false;
            warning( "Can't run parallel processes with GPU enabled. useGPU set to 'false'." )
        else
            opt.useGPU = oc.useGPU;
        end
    else
        opt.useGPU = oc.useGPU;
    end
else
    opt.useGPU = false;
    oc.useGPU = false;
    warning( "Setting useGPU to 'false'." )
end

if ( gpuDeviceCount < 1 )
    opt.useGPU = false;
    oc.useGPU = false;
    warning( "Setting useGPU to 'false'." )
end

if strcmpi( opt.structtype, "val" )
    opt.useGPU = false;
    oc.useGPU = false;
end


% fprintf( "\n----------------------\n" );
% fprintf( "assignGPUInformation\n\n" );
% fprintf( "opt.useGPU:\t%s\n", string( opt.useGPU ) );
% fprintf( "oc.useGPU:\t%s", string( oc.useGPU ) );
% fprintf( "\n----------------------\n" );

end