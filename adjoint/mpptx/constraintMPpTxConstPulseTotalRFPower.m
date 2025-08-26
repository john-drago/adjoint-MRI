function [ c, gradc ] = constraintMPpTxConstPulseTotalRFPower( pSc, opt )
% This function will calculate the RF power constraint value and (possibly
% the gradient) based on the design vector, pSc.

% Assume p is scaled to be between -1 and 1
pSc = pSc( : );

totalRFPower_constr = opt.totalRFPower_constr;

dutyCycle = opt.dutyCycle;
numTimePoints = opt.numTimePoints;
RFSlewIntPtNum = opt.RFSlewIntPtNum;
% numXYCoils = opt.numXYCoils;
% num_ORSP = opt.num_ORSP;
% num_MPSP = opt.num_MPSP;
tORSPIntPtNum = opt.tORSPIntPtNum;
tMPSPIntPtNum = opt.tMPSPIntPtNum;

slewScale = ( ( -0.5 + (1:RFSlewIntPtNum) ) / RFSlewIntPtNum );

brealORSP_idx = opt.breal_ORSP_idx;
brealORSP_sc = opt.scVec( brealORSP_idx );
brealORSP = brealORSP_sc .* pSc( brealORSP_idx );
brealORSP_pwr = brealORSP.^2;
c_brealORSP_pwr = tORSPIntPtNum * brealORSP_pwr;
brealORSP_slew = repmat( brealORSP, [ 1, RFSlewIntPtNum ] ) .* ( slewScale );
brealORSP_slew_pwr_sum = sum( brealORSP_slew.^2, 2 );
c_brealORSP_slew_pwr = 2 * sum( brealORSP_slew_pwr_sum, 2 );

brealMPSP_idx = opt.breal_MPSP_idx;
brealMPSP_sc = opt.scVec( brealMPSP_idx );
brealMPSP = brealMPSP_sc .* pSc( brealMPSP_idx );
brealMPSP_pwr = brealMPSP.^2;
c_brealMPSP_pwr = tMPSPIntPtNum * brealMPSP_pwr;
brealMPSP_slew = repmat( brealMPSP, [ 1, RFSlewIntPtNum ] ) .* ( slewScale );
brealMPSP_slew_pwr_sum = sum( brealMPSP_slew.^2, 2 );
c_brealMPSP_slew_pwr = 2 * sum( brealMPSP_slew_pwr_sum, 2 );

bimagORSP_idx = opt.bimag_ORSP_idx;
bimagORSP_sc = opt.scVec( bimagORSP_idx );
bimagORSP = bimagORSP_sc .* pSc( bimagORSP_idx );
bimagORSP_pwr = bimagORSP.^2;
c_bimagORSP_pwr = tORSPIntPtNum * bimagORSP_pwr;
bimagORSP_slew = repmat( bimagORSP, [ 1, RFSlewIntPtNum ] ) .* ( slewScale );
bimagORSP_slew_pwr_sum = sum( bimagORSP_slew.^2, 2 );
c_bimagORSP_slew_pwr = 2 * sum( bimagORSP_slew_pwr_sum, 2 );

bimagMPSP_idx = opt.bimag_MPSP_idx;
bimagMPSP_sc = opt.scVec( bimagMPSP_idx );
bimagMPSP = bimagMPSP_sc .* pSc( bimagMPSP_idx );
bimagMPSP_pwr = bimagMPSP.^2;
c_bimagMPSP_pwr = tMPSPIntPtNum * bimagMPSP_pwr;
bimagMPSP_slew = repmat( bimagMPSP, [ 1, RFSlewIntPtNum ] ) .* ( slewScale );
bimagMPSP_slew_pwr_sum = sum( bimagMPSP_slew.^2, 2 );
c_bimagMPSP_slew_pwr = 2 * sum( bimagMPSP_slew_pwr_sum, 2 );

c_unsc = dutyCycle/( 2 * opt.Z0 * numTimePoints ) *...
    ( sum(...
    c_brealORSP_pwr + c_brealORSP_slew_pwr + c_brealMPSP_pwr + c_brealMPSP_slew_pwr +...
    c_bimagORSP_pwr + c_bimagORSP_slew_pwr + c_bimagMPSP_pwr + c_bimagMPSP_slew_pwr ) )...
    - totalRFPower_constr;
c = c_unsc / totalRFPower_constr;

if nargout > 1
    
    % Initialize gradient vector
    gradc_unsc = zeros( opt.numVars, 1 );

    % determine sum of points in the slew ramp
    sumSlewRamp = sum( slewScale.^2 );

    % Assign gradients with respect to breal (or bx)
    grad_brealORSP_slew = 2 * sumSlewRamp * brealORSP;
    grad_brealORSP = 2 * brealORSP;
    gradc_unsc( brealORSP_idx, 1 ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( tORSPIntPtNum * grad_brealORSP + 2 * grad_brealORSP_slew ) .* brealORSP_sc;

    grad_brealMPSP_slew = 2 * sumSlewRamp * brealMPSP;
    grad_brealMPSP = 2 * brealMPSP;
    gradc_unsc( brealMPSP_idx, 1 ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( tMPSPIntPtNum * grad_brealMPSP + 2 * grad_brealMPSP_slew ) .* brealMPSP_sc;
    
    % Assign gradients with respect to bimag (or by)
    grad_bimagORSP_slew = 2 * sumSlewRamp * bimagORSP;
    grad_bimagORSP = 2 * bimagORSP;
    gradc_unsc( bimagORSP_idx, 1 ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( tORSPIntPtNum * grad_bimagORSP + 2 * grad_bimagORSP_slew ) .* bimagORSP_sc;

    grad_bimagMPSP_slew = 2 * sumSlewRamp * bimagMPSP;
    grad_bimagMPSP = 2 * bimagMPSP;
    gradc_unsc( bimagMPSP_idx, 1 ) = ( ( dutyCycle ) / ( 2 * opt.Z0 * numTimePoints ) ) *...
        ( tMPSPIntPtNum * grad_bimagMPSP + 2 * grad_bimagMPSP_slew ) .* bimagMPSP_sc;

    gradc = gradc_unsc / totalRFPower_constr;

end
end