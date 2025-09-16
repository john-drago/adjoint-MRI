function [ A_acc, b_acc ] = generateAccelMatrix( varName, opt, accelConstr )

% Create "positive" acceleration matrix
acc_i = zeros( opt.varAmts( opt.varNumDictionary( varName ) )*3, 1 );
acc_j = zeros( opt.varAmts( opt.varNumDictionary( varName ) )*3, 1 );
acc_v = zeros( opt.varAmts( opt.varNumDictionary( varName ) )*3, 1 );

acc_idxs = opt.varIdxs( opt.varNumDictionary( varName ) );
acc_idxs = acc_idxs{ 1 };

sc_vals = [...
    opt.scVec( acc_idxs( 1 ) ); opt.scVec( acc_idxs ); opt.scVec( acc_idxs( end ) ) ];

% -2 on the diagonal
acc_i( 1:opt.varAmts( opt.varNumDictionary( varName ) ) ) =...
    1 + ( 1:opt.varAmts( opt.varNumDictionary( varName ) ) ).';
acc_j( 1:opt.varAmts( opt.varNumDictionary( varName ) ) ) =...
    acc_idxs;
acc_v( 1:opt.varAmts( opt.varNumDictionary( varName ) ) ) =...
    -2 * ones( size( acc_idxs ) );

% 1 on first superdiagonal
acc_i( ( opt.varAmts( opt.varNumDictionary( varName ) ) + 1 ):( 2 * opt.varAmts( opt.varNumDictionary( varName ) ) ) ) =...
    ( 1 : ( opt.varAmts( opt.varNumDictionary( varName ) ) ) ).';
acc_j( ( opt.varAmts( opt.varNumDictionary( varName ) ) + 1 ):( 2 * opt.varAmts( opt.varNumDictionary( varName ) ) ) ) =...
    acc_idxs( 1:end );
acc_v( ( opt.varAmts( opt.varNumDictionary( varName ) ) + 1 ):( 2 * opt.varAmts( opt.varNumDictionary( varName ) ) ) ) =...
    1 * ones( ( length( acc_idxs ) ), 1 );

% 1 on first subdiagonal
acc_i( ( 2 * opt.varAmts( opt.varNumDictionary( varName ) ) + 1 ):( 3 * opt.varAmts( opt.varNumDictionary( varName ) ) ) ) =...
    2 + ( 1: ( opt.varAmts( opt.varNumDictionary( varName ) ) ) ).';
acc_j( ( 2 * opt.varAmts( opt.varNumDictionary( varName ) ) + 1 ):( 3 * opt.varAmts( opt.varNumDictionary( varName ) ) ) ) =...
    acc_idxs( 1:(end) );
acc_v( ( 2 * opt.varAmts( opt.varNumDictionary( varName ) ) + 1 ):( 3 * opt.varAmts( opt.varNumDictionary( varName ) ) ) ) =...
    1 * ones( ( length( acc_idxs )  ), 1  );

% Generate matrices
A_acc_pos = sparse( acc_i, acc_j, acc_v,...
    opt.varAmts( opt.varNumDictionary( varName ) )+2, opt.numVars );
% Generate "negative" difference matrix for breal_ORSP
A_acc_neg = -1 * A_acc_pos;

% Concatenate matrices
A_acc = [ A_acc_pos; A_acc_neg ];
b_acc = ( ( opt.dt^2 * accelConstr ) ./ [ sc_vals; sc_vals ] ) .*...
    ones( 2 * ( opt.varAmts( opt.varNumDictionary( varName ) ) + 2 ), 1 );
end