function [ partialfpartialp, varInfo, initSt ] = ...
    forwardModelGradientCheb_base( partialfpartialp, initSt, wv, opt, nn )
% This function will calculate the gradp_fArray at a given time interval.

%% Calculate partialfpartialB function first
% initSt = handlepartialfpartialBRotationMatrix( opt, initSt, wv, Mtt, Mttp1, tt );

%% Initialize varInfo struct
varInfo = initializeVarInfoForwardModelGradient( opt );

%% Get partialB partial parameter for the given time point 
% determine what phase of the pulse we are in

% get time
% ti = opt.tvec( nn );
% dti = opt.dtvec( nn );

%% RF
% Calculate contribution from breal and bimag
[ partialfpartialbrealcheb, partialfpartialbimagcheb ] =...
    calcpartialBpartialRFCheb( initSt.partialfpartialBx, initSt.partialfpartialBy, opt, wv, nn );

% Add contribution from breal
varInfo = addVarIdxs( opt.breal_idx, varInfo );
partialfpartialp = placepartialfpartialpCalculation(...
    partialfpartialp,...
    partialfpartialbrealcheb,...
    varInfo );

% Add contribution from bimag
varInfo = addVarIdxs( opt.bimag_idx, varInfo );
partialfpartialp = placepartialfpartialpCalculation(...
    partialfpartialp,...
    partialfpartialbimagcheb,...
    varInfo );

%% Grad
varInfo = addVarIdxs( opt.grad_idx, varInfo );
partialfpartialp = placepartialfpartialpCalculation(...
    partialfpartialp,...
    calcpartialBpartialGradCheb( initSt.partialfpartialBz, opt, wv, nn ),...
    varInfo );

%% Shim
if opt.numZCoils > 0
    varInfo = addVarIdxs( opt.shim_idx, varInfo );
    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        calcpartialBpartialShimCheb( initSt.partialfpartialBz, opt, wv, nn ),...
        varInfo );
end

end