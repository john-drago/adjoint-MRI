function VOPs = performVOPCompression_Lee( QLocalMats, epsG, ovMat )
% This function will perform VOP Compression for the QLocalMats that are
% passed in according to the Lee algorithm for a given overestimation
% factor epsilonG and an overestimation matrix, ovMat.
%
% Lee, MRM, 2012. doi: 10.1002/mrm.23140
%
% The implementation does NOT perform the iterative algorithm specified in
% step 3 of the paper (pg. 1569), but it instead performs an optimization
% using gradient descent that was mentioned by Kuehne 2017 ISMRM during
% presentation.
% https://cds.ismrm.org/protected/17MPresentations/videos/0478/index.html
%
% The full Kuehne idea has not been implemented however

%% Determine maximum eigenvalues for all of the local SAR matrices
numVox = size( QLocalMats, 3 );
maxEigQVox = zeros( numVox, 1 );
minEigQVox = zeros( numVox, 1 );

for nn = 1:numVox

    % Ensure that the diagonal elements are real (condition for PSD matrix)
    QLocalMats( :, :, nn ) = QLocalMats( :, :, nn )...
        - diag( diag( QLocalMats( :, :, nn ) ) ) + diag( real( diag( QLocalMats( :, :, nn ) ) ) );

    eigQVox = eig( QLocalMats( :, :, nn ) );

    maxEigQVox( nn ) = max( real( eigQVox ) );
    minEigQVox( nn ) = min( real( eigQVox ) );

end

%% Sort Eigenvalues in descending order
[ maxEigQVox, maxEigidx ] = sort( maxEigQVox, 'descend' );

%% Determine Overestimation Matrix
if isnumeric( ovMat )
    ovMatInit = ovMat;
elseif ischar( ovMat )
    if matches( ovMat, 'identity', 'IgnoreCase', true )
        ovMatInit = eye( size( QLocalMats, [ 1 2 ] ) );
    elseif matches( ovMat, 'maxlocalSAReig', 'IgnoreCase', true )
        ovMatInit = maxEigQVox(1) * eye( size( QLocalMats, [ 1 2 ] ) );
    end
end

% calculate the final overestimation matrix that depends on both the matrix
% that is passed in and epsG
ovMat = ovMatInit * epsG;

%% Make the set of VOPs and subset of Q: Qsub
% initialize arrays
Qsub_log = false( numVox, 1 );
Qstilltoeval_log = true( numVox, 1 );
PSDfeas_log = true( numVox, 1 );

% We already know that the first entry is the first VOP
numVOP = 1;
Qsub_log( 1 ) = true;
Qstilltoeval_log( 1 ) = false;
PSDfeas_log( 1 ) = false;

for mm = 2:numVox

    % Get current matrix under consideration
    Qcurrent = QLocalMats( :, :, maxEigidx( mm ) );
    Qstilltoeval_log( mm ) = false;

    % Get current Qsub
    Qsub = QLocalMats( :, :, maxEigidx( Qsub_log ) );

    % Determine if it is possible to bound this matrix with other matrices
    % assuming naive c_wv:
    c_wv_0 = ones( numVOP, 1 ) / numVOP;
    c_wv_0 = reshape( c_wv_0, [ 1, 1, length(c_wv_0) ] );
    LMI_0 = sum( ( Qsub ) .* c_wv_0, 3 ) + ovMat - Qcurrent;
    [ ~, cholflag_0 ] = chol( LMI_0 ); % use cholesky as this is a quick test for PSD
    
    % Determine if it is possible to bound this matrix with the previously
    % optimized c_wv
    if ( numVOP > 1 ) && ( cholflag_0 ~= 0 )
        c_wv_feas = reshape( c_wv_feas, [ 1, 1, length(c_wv_feas) ] );
        if ~converged
            % determines what the size of c_wv_feas should be
            LMI_prev = sum( ( Qsub( :, :, 1:(end-1) ) ) .* c_wv_feas, 3 ) + ovMat - Qcurrent;
        else
            % c_wv_feas will be one element larger in this case
            LMI_prev = sum( ( Qsub ) .* c_wv_feas, 3 ) + ovMat - Qcurrent;
        end
        [ ~, cholflag_prev ] = chol( LMI_prev ); % use cholesky as this is a quick test for PSD
    else
        cholflag_prev = 1;
    end

    % If we can't bound matrix, then proceed with attempt to find feasible point
    if ( cholflag_0 ~= 0 ) && ( cholflag_prev ~= 0 )
        
        % run the feasibility problem of Lee algorithm
        [ c_wv_feas, converged ] = feasibilityProblemLee(...
            squeeze( c_wv_0 ), Qcurrent, Qsub, ovMat );
        if ~converged
            numVOP = numVOP + 1;
            Qsub_log( mm ) = true; % mark this Qsub as a VOP
        end
    else % we can bound the matrix naively 
        PSDfeas_log( mm ) = false; % if we don't have to perform feasibility problem
    end

end

%% Assemble the VOPs
VOPs = QLocalMats( :, :, maxEigidx( Qsub_log ) ) + ovMat; % make sure to add the ovMat factor

end