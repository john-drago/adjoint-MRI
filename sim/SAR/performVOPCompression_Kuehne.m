function [ VOPs, VOPID ] = performVOPCompression_Kuehne( QLocalMats, epsG, ovMat )
% This function will perform VOP Compression for the QLocalMats that are
% passed in according to the Kuehne adaptation of the Lee algorithm for a
% given overestimation factor epsilonG and an overestimation matrix, ovMat.
%
% Kuehne, ISMRM 2017.
% https://cds.ismrm.org/protected/17MPresentations/videos/0478/index.html
%
% The implementation does NOT perform the iterative algorithm specified in
% step 3 of the Lee paper (doi: 10.1002/mrm.23140), but it instead performs
% an optimization using gradient descent that was mentioned by Kuehne 2017
% ISMRM during presentation.
%

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
PSDopt_log = false( numVox, 1 );

% We already know that the first entry is the first VOP
numVOP = 1;
Qsub_log( 1 ) = true;
Qstilltoeval_log( 1 ) = false;
numQtoEval = sum( Qstilltoeval_log );

mm = 2;

while ( numQtoEval > 0 )

    % Get current matrix under consideration
    Qcurrent_mm = QLocalMats( :, :, maxEigidx( mm ) );
    Qstilltoeval_log( mm ) = false;

    % Get current Qsub
    Qsub = QLocalMats( :, :, maxEigidx( Qsub_log ) );

    % Determine if it is possible to bound this matrix with other matrices
    % assuming naive c_wv:
    c_wv_0 = ones( numVOP, 1 ) / numVOP;
    LMI_0 = sum( ( Qsub ) .* reshape( c_wv_0, [ 1, 1, length(c_wv_0) ] ), 3 ) + ovMat - Qcurrent_mm;
    [ ~, cholflag_0 ] = chol( LMI_0 ); % use cholesky as this is a quick test for PSD

    % if can't bound matrix naively, run optimization of Kuehne
    if ( cholflag_0 ~= 0 )

        [ c_wv_opt, out ] = optProblemKuehne(...
            squeeze( c_wv_0 ), Qcurrent_mm, Qsub, ovMat );
        PSDopt_log( mm ) = true;

        if out > 0
            numVOP = numVOP + 1;
            Qsub_log( mm ) = true; % mark this Qsub as a VOP

            mm = find( Qstilltoeval_log, 1, 'first' );
            numQtoEval = numQtoEval - 1;
        else

            LMI_pre = 1/(sum(c_wv_opt)) * sum( ( Qsub ) .* reshape( c_wv_opt, [ 1, 1, length(c_wv_opt) ] ), 3 ) + ovMat;
            
            % could try to vectorize a method to determine if diagonal of
            % LMI_pre for all future matrices is strictly positive. For those
            % that are not strictly positive, we know it won't be PSD.

            % determine voxels left
            Qstilltoeval_log_rest_mm = Qstilltoeval_log( (mm+1):numVox );
            Qstilltoeval_idx_rest_mm = (mm+1):numVox;
            Qstilltoeval_idx_mm = Qstilltoeval_idx_rest_mm( Qstilltoeval_log_rest_mm );
            
            firstUndominated = true;
            for nn = 1:length( Qstilltoeval_idx_mm )
                Qcurrent_nn = QLocalMats( :, :, maxEigidx( Qstilltoeval_idx_mm( nn ) ) );
                LMI_nn = LMI_pre - Qcurrent_nn;

                [ ~, cholflag_nn ] = chol( LMI_nn ); % use cholesky as this is a quick test for PSD
                if cholflag_nn == 0
                    Qstilltoeval_log( Qstilltoeval_idx_mm( nn ) ) = false;
                else
                    if firstUndominated
                        mm = Qstilltoeval_idx_mm( nn );
                        firstUndominated = false;
                    end
                end
            end
            
            numQtoEval = sum( Qstilltoeval_log );
        end
    else
        mm = find( Qstilltoeval_log, 1, 'first' );
        numQtoEval = numQtoEval - 1;
    end
end

%% Assemble the VOPs
VOPs = QLocalMats( :, :, maxEigidx( Qsub_log ) ) + ovMat; % make sure to add the ovMat factor
VOPID = maxEigidx( Qsub_log );

end