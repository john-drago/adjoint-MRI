function [oc, opt] = assignGPUInformation(oc, opt)
% Decide whether to use GPU based on user request, optimizer parallelism,
% device availability, and struct type. Keeps oc.useGPU and opt.useGPU in sync.

    % Helpers
    function tf = getUseParallel(st, field)
        tf = false;
        if isfield(st, field) && isstruct(st.(field)) && isfield(st.(field), 'UseParallel')
            tf = logical(st.(field).UseParallel);
        end
    end

    function [avail,count] = gpuAvail()
        avail = false; count = 0;
        % Prefer robust check; handle missing toolbox gracefully
        if exist('gpuDeviceCount','file') == 2
            try
                count = gpuDeviceCount;
                avail = count > 0;
            catch
                avail = false; count = 0;
            end
        elseif exist('parallel.gpu.GPUDevice','class') == 8
            try
                avail = parallel.gpu.GPUDevice.isAvailable;
                count = double(avail); % unknown exact count here
            catch
                avail = false; count = 0;
            end
        end
    end

    % Initialize requested flag
    if ~isfield(oc, 'useGPU') || isempty(oc.useGPU)
        oc.useGPU = false;
        % (No warning here; being default is fine.)
    end

    % Optimizer settings
    fminUseParallel = getUseParallel(oc, 'fminopt');
    gaUseParallel   = getUseParallel(oc, 'gaopt');  % <-- fixed

    % Normalize optType for comparisons
    optType = string(getfield(oc, 'optType')); %#ok<GFLD>
    isGaProxy = strcmpi(optType, "ga-proxy");

    % Start with requested value
    opt.useGPU = logical(oc.useGPU);

    % Disallow GPU if structtype == "val"
    if isfield(opt, 'structtype') && strcmpi(string(opt.structtype), "val")
        if opt.useGPU
            warning("assignGPUInformation:StructTypeVal", ...
                "useGPU disabled because opt.structtype is 'val'.");
        end
        opt.useGPU = false;
    end

    % Disallow GPU + Parallel clash (except ga-proxy passthrough)
    if ~isGaProxy
        if (fminUseParallel || gaUseParallel) && oc.useGPU
            warning("assignGPUInformation:ParallelClash", ...
                "Parallel optimization is enabled; disabling GPU (cannot use GPU and parallel workers together).");
            opt.useGPU = false;
        end
    end

    % Disable if no GPU devices available
    [gpuOK, nGPU] = gpuAvail();
    if opt.useGPU && ~gpuOK
        warning("assignGPUInformation:NoGPU", ...
            "No GPU device detected (gpuDeviceCount=%d). Disabling GPU.", nGPU);
        opt.useGPU = false;
    end

    % Keep oc and opt consistent 
    oc.useGPU = opt.useGPU;

    % % Uncomment for debug prints
    % fprintf("\n----------------------\nassignGPUInformation\n");
    % fprintf("opt.useGPU:\t%s\n", string(opt.useGPU));
    % fprintf("oc.useGPU:\t%s\n", string(oc.useGPU));
    % fprintf("----------------------\n");
end
