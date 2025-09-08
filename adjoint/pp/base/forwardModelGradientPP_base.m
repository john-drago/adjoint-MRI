function [ partialfpartialp, varInfo, initSt ] = ...
    forwardModelGradientPP_base( partialfpartialp, initSt, wv, opt, nn )
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
% determine RF period
PPRFidx = all( ( nn >= opt.PPIdxs_RF( :, 1 ) ) & ( nn <= opt.PPIdxs_RF( :, 2 ) ), 2 );

if any( PPRFidx )

    % Create shapeFnValsTimePoints_RF_rep
    shapeFnValsTimePoints_RF_rep = repmat( ...
        wv.varsToTimepoints_RF( nn, : ), [ opt.numPos, 1 ] );

    PPRFVaridx = opt.PPVarIdxsCell_RF{ PPRFidx };

    % Calculate contribution from breal and bimag
    [ partialfpartialbrealshape, partialfpartialbimagshape ] =...
        calcpartialBpartialRFPP( initSt.partialfpartialBx, initSt.partialfpartialBy,...
        shapeFnValsTimePoints_RF_rep( :, PPRFVaridx ), opt, wv  );

    % Add contribution from breal
    varInfo = addVarIdxs( opt.breal_idx( opt.PPVarIdxsCoilCell_RF{ PPRFidx } ), varInfo );
    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        partialfpartialbrealshape,...
        varInfo );

    % Add contribution from bimag
    varInfo = addVarIdxs( opt.bimag_idx( opt.PPVarIdxsCoilCell_RF{ PPRFidx } ), varInfo );
    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        partialfpartialbimagshape,...
        varInfo );
end

%% Grad
% determine Grad period
PPgradidx = all( ( nn >= opt.PPIdxs_grad( :, 1 ) ) & ( nn <= opt.PPIdxs_grad( :, 2 ) ), 2 );

if any( PPgradidx )

    % Create shapeFnValsTimePoints_grad_rep
    shapeFnValsTimePoints_grad_rep = repmat( ...
        wv.varsToTimepoints_grad( nn, : ), [ opt.numPos, 1 ] );

    PPgradVaridx = opt.PPVarIdxsCell_grad{ PPgradidx };

    varInfo = addVarIdxs( opt.grad_idx( opt.PPVarIdxsCoilCell_grad{ PPgradidx } ), varInfo );
    partialfpartialp = placepartialfpartialpCalculation(...
        partialfpartialp,...
        calcpartialBpartialGradPP( initSt.partialfpartialBz,...
        shapeFnValsTimePoints_grad_rep( :, PPgradVaridx ),...
        opt, wv  ),...
        varInfo );

end

%% Shim
if opt.numZCoils > 0
    % determine Shim period
    PPshimidx = all( ( nn >= opt.PPIdxs_shim( :, 1 ) ) & ( nn <= opt.PPIdxs_shim( :, 2 ) ), 2 );

    if any( PPshimidx )
        PPshimVaridx = opt.PPVarIdxsCell_shim{ PPshimidx };

        % Create shapeFnValsTimePoints_shim_rep
        shapeFnValsTimePoints_shim_rep = repmat( ...
            wv.varsToTimepoints_shim( nn, : ), [ opt.numPos, 1 ] );

        varInfo = addVarIdxs( opt.shim_idx( opt.PPVarIdxsCoilCell_shim{ PPshimidx } ), varInfo );
        partialfpartialp = placepartialfpartialpCalculation(...
            partialfpartialp,...
            calcpartialBpartialShimPP( initSt.partialfpartialBz,...
            shapeFnValsTimePoints_shim_rep( :, PPshimVaridx ),...
            opt, wv  ),...
            varInfo );
    end
end

end