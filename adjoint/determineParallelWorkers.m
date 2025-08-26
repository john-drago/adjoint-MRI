function oc = determineParallelWorkers( numWorkers, oc )

if isfield( oc, 'fminopt' )
    fminoptUseParallel = oc.fminopt.UseParallel;
else
    fminoptUseParallel = false;
end

if isfield( oc, 'gaopt' )
    gaoptUseParallel = oc.gaopt.UseParallel;
else
    gaoptUseParallel = false;
end

pp = gcp('nocreate');
if fminoptUseParallel || gaoptUseParallel
    numcores = feature('numcores');
    if isempty(pp)
        pp = parpool( min( [ numWorkers, numcores ] ) );
    else
        if pp.NumWorkers ~= numWorkers
            delete( pp );
            pp = parpool( min( [ numWorkers, numcores ] ) );
        end
    end
    oc.numWorkers = pp.NumWorkers;
else
    delete( pp );
    oc.numWorkers = 1;
end

end