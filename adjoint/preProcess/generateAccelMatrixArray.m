function [ A_acc, b_acc ] = generateAccelMatrixArray( varName, numVarType, opt, accelConstr )

varAmts = opt.varAmts( opt.varNumDictionary( varName ) );
varAmts_nn = varAmts / numVarType;

acc_idxs = opt.varIdxs( opt.varNumDictionary( varName ) );
acc_idxs = acc_idxs{ 1 };

numRows = numVarType * ( varAmts_nn + 2 );

% Create "positive" acceleration matrix
acc_i = zeros( numVarType * (varAmts_nn * 3), 1 );
acc_j = zeros( numVarType * (varAmts_nn * 3), 1 );
acc_v = zeros( numVarType * (varAmts_nn * 3), 1 );
sc_vals = zeros( numRows, 1 );

for nn = 1:numVarType
    
    acc_i_nn = zeros( varAmts_nn * 3, 1 );
    acc_j_nn = zeros( varAmts_nn * 3, 1 );
    acc_v_nn = zeros( varAmts_nn * 3, 1 );

    acc_idxs_nn = acc_idxs( ( (nn - 1)*varAmts_nn + (1 : varAmts_nn) ) );

    sc_vals_nn = [ opt.scVec( acc_idxs_nn( 1 ) ); opt.scVec( acc_idxs_nn ) ; opt.scVec( acc_idxs_nn( end ) ) ];

    % -2 on the diagonal
    acc_i_nn( 1:varAmts_nn ) =...
        ( (nn - 1)*(varAmts_nn+2) + ( 1 + ( 1:varAmts_nn ) ) ).';
    acc_j_nn( 1:varAmts_nn ) =...
        acc_idxs_nn;
    acc_v_nn( 1:varAmts_nn ) =...
        -2 * ones( varAmts_nn,1 );

    % 1 on first superdiagonal
    acc_i_nn( ( varAmts_nn + 1 ):( 2 * varAmts_nn ) ) =...
        ( (nn - 1)*(varAmts_nn+2) + ( 1 : ( varAmts_nn ) ) ).';
    acc_j_nn( ( varAmts_nn + 1 ):( 2 * varAmts_nn ) ) =...
        acc_idxs_nn( 1:end );
    acc_v_nn( ( varAmts_nn + 1 ):( 2 * varAmts_nn ) ) =...
        1 * ones( (varAmts_nn), 1 );

    % 1 on first subdiagonal
    acc_i_nn( ( 2 * varAmts_nn + 1 ):( 3 * varAmts_nn ) ) =...
        ( (nn - 1)*(varAmts_nn+2) + 2 + ( 1: varAmts_nn ) ).';
    acc_j_nn( ( 2 * varAmts_nn + 1 ):( 3 * varAmts_nn ) ) =...
        acc_idxs_nn( 1:(end) );
    acc_v_nn( ( 2 * varAmts_nn + 1 ):( 3 * varAmts_nn ) ) =...
        1 * ones( (varAmts_nn), 1  );

    % Assign position
    acc_i( (nn - 1)*(3*varAmts_nn) + (1:(3*varAmts_nn)).' ) = acc_i_nn;
    acc_j( (nn - 1)*(3*varAmts_nn) + (1:(3*varAmts_nn)).' ) = acc_j_nn;
    acc_v( (nn - 1)*(3*varAmts_nn) + (1:(3*varAmts_nn)).' ) = acc_v_nn;

    sc_vals( (nn - 1)*( varAmts_nn+2 ) + ( 1:( varAmts_nn + 2 ) ).' ) =...
        sc_vals_nn;

end

% Generate matrices
A_acc_pos = sparse( acc_i, acc_j, acc_v,...
    numRows, opt.numVars );
% Generate "negative" difference matrix for breal_ORSP
A_acc_neg = -1 * A_acc_pos;

% Concatenate matrices
A_acc = [ A_acc_pos; A_acc_neg ];
b_acc = ( ( opt.dt^2 * accelConstr ) ./ [ sc_vals; sc_vals ] ) .*...
    ones( 2 * numRows, 1 );
end