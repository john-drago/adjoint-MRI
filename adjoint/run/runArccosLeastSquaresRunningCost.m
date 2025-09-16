function [ sum_g, g, gradx_g, gradp_g ] = runArccosLeastSquaresRunningCost( Marray, ~, opt )
% Vectorized running cost with arccos least-squares
%
% Handles:
%   - GPU mode: Marray is [ numPos, 3, numTimePoints+1 ]
%   - CPU/cell mode: Marray is { numTimePoints x 1 }, each cell [ numPos, 3 ]
%
% Outputs:
%   sum_g   : scalar, (sum over time, mean over subjects) * opt.runningCostWeightingCoefficient
%   g       : [ 1, numTimePoints+1 ] mean across subjects, per-time NRMSE (unscaled)
%   gradx_g : [ numPos, 3, numTimePoints+1 ] gradient w.r.t. state (nonzero only in z)
%   gradp_g : [ opt.numVars, 1 ] zeros (no direct param dependence)

numSubj = opt.numSubj;
numPos = opt.numPos;
numTimePoints = opt.numTimePoints;

% Ensure everyting is a single [ numPos, 3, numTimePoints+1 ] array without loops

if ~opt.useGPU && iscell( Marray )
    % CPU: transition cell array to numeric array
    likeT = Marray{ 1 };
    Marray = cat( 3, Marray{ : } ); % Stack cells into a 3-D array: [ numPos, 3 , numTimePoints ]
else
    likeT = Marray( :, :, 1 );
end

% Subject mapping and grouping matrix
% sIdx(i) = subject id for state i
sIdx = repelem( transpose( 1:numSubj ), opt.numPosSubj(:) );
sIdx = sIdx( : );
subjMat = sparse( 1:numPos, sIdx, 1, numPos, numSubj ); % [ numPos x numSubjects ] (column s has 1s for idxs of subject s)
if opt.useGPU
    subjMat = gpuArray( subjMat );
end

% Bound Mz into (-1,1) to keep acos real
Mz = reshape( Marray( :, 3, : ), [ numPos, (numTimePoints+1) ] );
oneLike = cast( 1, 'like', likeT );
epsLike = max( eps(oneLike), cast( 1e-12, 'like', likeT ) );
Mz = min( max( Mz, -oneLike + epsLike ), oneLike - epsLike );

% Target FA (time-invariant)
Mztarg = cast( opt.Mtarg( :, 3 ), 'like', likeT );
Mztarg = min( max(Mztarg, -oneLike + epsLike), oneLike - epsLike );
targFA = real( acos(Mztarg) );      
targFA2 = targFA.^2;
if opt.useGPU
    targFA2 = gpuArray( targFA2 );
    targFA = gpuArray( targFA );
end

% Subject-wise ||targFA||_2
targFA2_sum_per_subj = transpose( subjMat ) * targFA2;
targFAnorm_per_subj = sqrt( targFA2_sum_per_subj );
mask_targFA = (targFAnorm_per_subj > 0);

% Errors across time
FA = real( acos( Mz ) );
errFA = FA - targFA;

% Subject-wise ||errFA||_2 at each time
errFAnorm_over_time_per_subj = sqrt( transpose( subjMat ) * ( errFA.^2 ) );
mask_err = (errFAnorm_over_time_per_subj > 0);

% Per-subject/per-time NRMSE (unscaled)
den_tFA = targFAnorm_per_subj + ~mask_targFA; % avoid 0-div; will be masked to 0
g_over_time_per_subj = errFAnorm_over_time_per_subj ./ den_tFA;
g_over_time_per_subj(~mask_targFA, :) = 0; % zero contribution when ‖targFA‖_2 = 0

% Scalar running cost: mean over subjects, sum over time, then scale
g = mean( g_over_time_per_subj, 1 );
sum_g = opt.runningCostWeightingCoefficient * sum( g );

if nargout > 2
    % Gradient w.r.t magnetization (nonzero only in z)
    denom_sqrt = sqrt( 1 - Mz.^2 );
    targFAnorm_per_pos = targFAnorm_per_subj( sIdx );
    mask_targFA_per_pos = mask_targFA( sIdx );

    errFAnorm_over_time_per_pos = errFAnorm_over_time_per_subj( sIdx, : );
    mask_err_over_time_per_pos = mask_err( sIdx, : );

    gradx_g_per_pos = - errFA ./ ( targFAnorm_per_pos .* errFAnorm_over_time_per_pos .* denom_sqrt );
    gradx_g_per_pos( ~mask_targFA_per_pos, : ) = 0;
    gradx_g_per_pos( ~mask_err_over_time_per_pos ) = 0;
    gradx_g_per_pos( ~isfinite( gradx_g_per_pos ) ) = 0;

    gradx_g_per_pos_sc = ( opt.runningCostWeightingCoefficient / numSubj) * gradx_g_per_pos;

    gradx_g = zeros( numPos, 3, numTimePoints+1, 'like', likeT );
    gradx_g( :, 3, : ) = reshape( gradx_g_per_pos_sc, [ numPos, 1 , numTimePoints+1 ] );

    if ~opt.useGPU
        % Current gradx_g is [ numPos, 3, numTimePoints+1 ] numeric (CPU). Split along time without loops.
        tmp = num2cell( gradx_g, [1 2] );   % -> 1 x 1 x (numTimePoints+1) cell, each cell is [ numPos x 3 ]
        tmp = squeeze( tmp );               % -> 1 x (numTimePoints+1)
        gradx_g = reshape( tmp, [numTimePoints+1, 1] );   % -> (numTimePoints+1) x 1, each cell [ numPos x 3 ]
    end

    % No param dependence
    % gradp_g = zeros( opt.numVars, ( opt.numTimePoints + 1 ) , 'like', likeT );
    gradp_g = [];
end

end
