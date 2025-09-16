function opt = cleanUpAdjointStructs( opt )
% This file will process the opt struct to delete parameters that don't
% need to be deleted before the opt struct is saved.

if isfield( opt, 'bzsens' )
    opt = rmfield( opt, 'bzsens' );
end
if isfield( opt, 'b1p' )
    opt = rmfield( opt, 'b1p' );
end
if isfield( opt, 'b1preal' )
    opt = rmfield( opt, 'b1preal' );
end
if isfield( opt, 'b1pimag' )
    opt = rmfield( opt, 'b1pimag' );
end
if isfield( opt, 'db0' )
    opt = rmfield( opt, 'db0' );
end
if isfield( opt, 'pos' )
    opt = rmfield( opt, 'pos' );
end
if isfield( opt, 'X' )
    opt = rmfield( opt, 'X' );
end
if isfield( opt, 'Y' )
    opt = rmfield( opt, 'Y' );
end
if isfield( opt, 'Z' )
    opt = rmfield( opt, 'Z' );
end
if isfield( opt, 'M0' )
    opt = rmfield( opt, 'M0' );
end
if isfield( opt, 'opt_roi_array' )
    opt = rmfield( opt, 'opt_roi_array' );
end
if isfield( opt, 'M_targ_array' )
    opt = rmfield( opt, 'M_targ_array' );
end
if isfield( opt, 'FA_targ_array' )
    opt = rmfield( opt, 'FA_targ_array' );
end
if isfield( opt, 'M_opt_array' )
    opt = rmfield( opt, 'M_opt_array' );
end
if isfield( opt, 'FA_opt_array' )
    opt = rmfield( opt, 'FA_opt_array' );
end
if isfield( opt, 'M_BCHP_array' )
    opt = rmfield( opt, 'M_BCHP_array' );
end
if isfield( opt, 'FA_BCHP_array' )
    opt = rmfield( opt, 'FA_BCHP_array' );
end
if isfield( opt, 'Mtarg' )
    opt = rmfield( opt, 'Mtarg' );
end
if isfield( opt, 'FAtarg' )
    opt = rmfield( opt, 'FAtarg' );
end
if isfield( opt, 'FABCHP' )
    opt = rmfield( opt, 'FABCHP' );
end
if isfield( opt, 'MBCHP' )
    opt = rmfield( opt, 'MBCHP' );
end
if isfield( opt, 'MBCHP' )
    opt = rmfield( opt, 'MBCHP' );
end
if isfield( opt, 'MBCHP' )
    opt = rmfield( opt, 'MBCHP' );
end
if isfield( opt, 'MThru' )
    opt = rmfield( opt, 'MThru' );
end
if isfield( opt, 'MtargThru' )
    opt = rmfield( opt, 'MtargThru' );
end
if isfield( opt, 'FAThru' )
    opt = rmfield( opt, 'FAThru' );
end
if isfield( opt, 'FAtargThru' )
    opt = rmfield( opt, 'FAtargThru' );
end
if isfield( opt, 'zThru' )
    opt = rmfield( opt, 'zThru' );
end
if isfield( opt, 'MInPlane' )
    opt = rmfield( opt, 'MInPlane' );
end
if isfield( opt, 'FAInPlane' )
    opt = rmfield( opt, 'FAInPlane' );
end
if isfield( opt, 'MvecInPlane' )
    opt = rmfield( opt, 'MvecInPlane' );
end
if isfield( opt, 'FAvecInPlane' )
    opt = rmfield( opt, 'FAvecInPlane' );
end
if isfield( opt, 'zInPlane' )
    opt = rmfield( opt, 'zInPlane' );
end
if isfield( opt, 'optroiInPlane' )
    opt = rmfield( opt, 'optroiInPlane' );
end
if isfield( opt, 'xySamplePtzVals' )
    opt = rmfield( opt, 'xySamplePtzVals' );
end
if isfield( opt, 'xySamplePtFAVals' )
    opt = rmfield( opt, 'xySamplePtFAVals' );
end
if isfield( opt, 'xySamplePtFAtargVals' )
    opt = rmfield( opt, 'xySamplePtFAtargVals' );
end
if isfield( opt, 'xySamplePtPhVals' )
    opt = rmfield( opt, 'xySamplePtPhVals' );
end
if isfield( opt, 'xySamplePtMxVals' )
    opt = rmfield( opt, 'xySamplePtMxVals' );
end
if isfield( opt, 'xySamplePtMyVals' )
    opt = rmfield( opt, 'xySamplePtMyVals' );
end
if isfield( opt, 'subjSlicesGloIdxs' )
    opt = rmfield( opt, 'subjSlicesGloIdxs' );
end
if isfield( opt, 'subjNSlicesGloIdxs' )
    opt = rmfield( opt, 'subjNSlicesGloIdxs' );
end
if isfield( opt, 'A' )
    opt = rmfield( opt, 'A' );
end
if isfield( opt, 'b' )
    opt = rmfield( opt, 'b' );
end
if isfield( opt, 'Aeq' )
    opt = rmfield( opt, 'Aeq' );
end
if isfield( opt, 'beq' )
    opt = rmfield( opt, 'beq' );
end
if isfield( opt, 'L' )
    opt = rmfield( opt, 'L' );
end
if isfield( opt, 'costFnOpt' )
    opt = rmfield( opt, 'costFnOpt' );
end
if isfield( opt, 'uniqueXYPos' )
    opt = rmfield( opt, 'uniqueXYPos' );
end
if isfield( opt, 'XYPosToUniqueXYPos' )
    opt = rmfield( opt, 'XYPosToUniqueXYPos' );
end
if isfield( opt, 'uniqueXYPosToXYPos' )
    opt = rmfield( opt, 'uniqueXYPosToXYPos' );
end
if isfield( opt, 'slicesMIdx' )
    opt = rmfield( opt, 'slicesMIdx' );
end
if isfield( opt, 'allSlicesMIdx' )
    opt = rmfield( opt, 'allSlicesMIdx' );
end
if isfield( opt, 'sliceIdxs' )
    opt = rmfield( opt, 'sliceIdxs' );
end
if isfield( opt, 'uniqueIdxs' )
    opt = rmfield( opt, 'uniqueIdxs' );
end
if isfield( opt, 'subjSlicesLocIdxs' )
    opt = rmfield( opt, 'subjSlicesLocIdxs' );
end
if isfield( opt, 'subjNSlicesLocIdxs' )
    opt = rmfield( opt, 'subjNSlicesLocIdxs' );
end
if isfield( opt, 'VOPs' )
    opt = rmfield( opt, 'VOPs' );
end
if isfield( opt, 'VOPs_E10m' )
    opt = rmfield( opt, 'VOPs_E10m' );
end
if isfield( opt, 'VOPs_H10m' )
    opt = rmfield( opt, 'VOPs_H10m' );
end
if isfield( opt, 'QGlobal' )
    opt = rmfield( opt, 'QGlobal' );
end
if isfield( opt, 'gpu_M0' )
    opt = rmfield( opt, 'gpu_M0' );
end
if isfield( opt, 'gpu_b1preal' )
    opt = rmfield( opt, 'gpu_b1preal' );
end
if isfield( opt, 'gpu_b1pimag' )
    opt = rmfield( opt, 'gpu_b1pimag' );
end
if isfield( opt, 'gpu_b1p' )
    opt = rmfield( opt, 'gpu_b1p' );
end
if isfield( opt, 'gpu_pos' )
    opt = rmfield( opt, 'gpu_pos' );
end
if isfield( opt, 'gpu_bzsens' )
    opt = rmfield( opt, 'gpu_bzsens' );
end
if isfield( opt, 'gpu_db0' )
    opt = rmfield( opt, 'gpu_db0' );
end
if isfield( opt, 'gpu_numZCoils' )
    opt = rmfield( opt, 'gpu_numZCoils' );
end
if isfield( opt, 'gpu_numXYCoils' )
    opt = rmfield( opt, 'gpu_numXYCoils' );
end
if isfield( opt, 'gpu_numPos' )
    opt = rmfield( opt, 'gpu_numPos' );
end
if isfield( opt, 'gpu_gyro' )
    opt = rmfield( opt, 'gpu_gyro' );
end
if isfield( opt, 'gpu_Tn' )
    opt = rmfield( opt, 'gpu_Tn' );
end
if isfield( opt, 'gpu_Tn_rep' )
    opt = rmfield( opt, 'gpu_Tn_rep' );
end
% if isfield( opt, 'Tn' )
%     opt = rmfield( opt, 'Tn' );
% end
if isfield( opt, 'Tn_rep' )
    opt = rmfield( opt, 'Tn_rep' );
end
if isfield( opt, 'gpu_FBRF' )
    opt = rmfield( opt, 'gpu_FBRF' );
end
if isfield( opt, 'gpu_FBgrad' )
    opt = rmfield( opt, 'gpu_FBgrad' );
end
if isfield( opt, 'gpu_FBshim' )
    opt = rmfield( opt, 'gpu_FBshim' );
end
if isfield( opt, 'costFnOpt' )
    opt = rmfield( opt, 'costFnOpt' );
end
if isfield( opt, 'nonlcon' )
    opt = rmfield( opt, 'nonlcon' );
end
% if isfield( opt, 'varsToTimepoints_RF' )
%     opt = rmfield( opt, 'varsToTimepoints_RF' );
% end
if isfield( opt, 'varsToChebByPeriods_RF' )
    opt = rmfield( opt, 'varsToChebByPeriods_RF' );
end
if isfield( opt, 'shapeFnChebCoeffs_RF' )
    opt = rmfield( opt, 'shapeFnChebCoeffs_RF' );
end
if isfield( opt, 'shapeFnValsTimePoints_RF_rep' )
    opt = rmfield( opt, 'shapeFnValsTimePoints_RF_rep' );
end
if isfield( opt, 'gpu_shapeFnValsTimePoints_RF_rep' )
    opt = rmfield( opt, 'gpu_shapeFnValsTimePoints_RF_rep' );
end
% if isfield( opt, 'varsToTimepoints_grad' )
%     opt = rmfield( opt, 'varsToTimepoints_grad' );
% end
if isfield( opt, 'varsToChebByPeriods_grad' )
    opt = rmfield( opt, 'varsToChebByPeriods_grad' );
end
if isfield( opt, 'shapeFnChebCoeffs_grad' )
    opt = rmfield( opt, 'shapeFnChebCoeffs_grad' );
end
if isfield( opt, 'shapeFnValsTimePoints_grad_rep' )
    opt = rmfield( opt, 'shapeFnValsTimePoints_grad_rep' );
end
if isfield( opt, 'gpu_shapeFnValsTimePoints_grad_rep' )
    opt = rmfield( opt, 'gpu_shapeFnValsTimePoints_grad_rep' );
end
% if isfield( opt, 'varsToTimepoints_shim' )
%     opt = rmfield( opt, 'varsToTimepoints_shim' );
% end
if isfield( opt, 'varsToChebByPeriods_shim' )
    opt = rmfield( opt, 'varsToChebByPeriods_shim' );
end
if isfield( opt, 'shapeFnChebCoeffs_shim' )
    opt = rmfield( opt, 'shapeFnChebCoeffs_shim' );
end
if isfield( opt, 'shapeFnValsTimePoints_shim_rep' )
    opt = rmfield( opt, 'shapeFnValsTimePoints_shim_rep' );
end
if isfield( opt, 'gpu_shapeFnValsTimePoints_shim_rep' )
    opt = rmfield( opt, 'gpu_shapeFnValsTimePoints_shim_rep' );
end
if isfield( opt, 'ipopt' )
    if isfield( opt.ipopt, 'funcs' )
        opt.ipopt = rmfield( opt.ipopt, 'funcs' );
    end
    if isfield( opt.ipopt, 'options' )
        if isfield( opt.ipopt.options, 'cl' )
            opt.ipopt.options = rmfield( opt.ipopt.options, 'cl' );
        end
        if isfield( opt.ipopt.options, 'cu' )
            opt.ipopt.options = rmfield( opt.ipopt.options, 'cu' );
        end
        if isfield( opt.ipopt.options, 'lb' )
            opt.ipopt.options = rmfield( opt.ipopt.options, 'lb' );
        end
        if isfield( opt.ipopt.options, 'ub' )
            opt.ipopt.options = rmfield( opt.ipopt.options, 'ub' );
        end
    end
end
if isfield( opt, 'nlopt' )
    if isfield( opt.nlopt, 'min_objective' )
        opt.nlopt = rmfield( opt.nlopt, 'min_objective' );
    end
    if isfield( opt.nlopt, 'mfc' )
        opt.nlopt = rmfield( opt.nlopt, 'mfc' );
    end
    if isfield( opt.nlopt, 'mfc_tol' )
        opt.nlopt = rmfield( opt.nlopt, 'mfc_tol' );
    end
    if isfield( opt.nlopt, 'mfc_count' )
        opt.nlopt = rmfield( opt.nlopt, 'mfc_count' );
    end
    if isfield( opt.nlopt, 'mh' )
        opt.nlopt = rmfield( opt.nlopt, 'mh' );
    end
    if isfield( opt.nlopt, 'mh_tol' )
        opt.nlopt = rmfield( opt.nlopt, 'mh_tol' );
    end
    if isfield( opt.nlopt, 'mh_count' )
        opt.nlopt = rmfield( opt.nlopt, 'mh_count' );
    end
end

if isfield( opt, 'nonlconSNOPT' )
    opt = rmfield( opt, 'nonlconSNOPT' );
end

if isfield( opt, 'wv' )
    if isfield( opt.wv, 'Tn_rep' )
        opt.wv = rmfield( opt.wv, 'Tn_rep' );
    end
    if isfield( opt.wv, 'shapeFnValsTimePoints_RF_rep' )
        opt.wv = rmfield( opt.wv, 'shapeFnValsTimePoints_RF_rep' );
    end
    if isfield( opt.wv, 'shapeFnValsTimePoints_grad_rep' )
        opt.wv = rmfield( opt.wv, 'shapeFnValsTimePoints_grad_rep' );
    end
    if isfield( opt.wv, 'shapeFnValsTimePoints_shim_rep' )
        opt.wv = rmfield( opt.wv, 'shapeFnValsTimePoints_shim_rep' );
    end
end

end