function [ A_slew, b_slew ] = generateSlewRateMatrix( varName, opt, slewConstr )
% Create "positive" difference matrix
% ones on the diagonal
slew_i = zeros( opt.varAmts( opt.varNumDictionary( varName ) )*2, 1 );
slew_j = zeros( opt.varAmts( opt.varNumDictionary( varName ) )*2, 1 );
slew_v = zeros( opt.varAmts( opt.varNumDictionary( varName ) )*2, 1 );

slew_idxs = opt.varIdxs( opt.varNumDictionary( varName ) );
slew_idxs = slew_idxs{ 1 };

sc_vals = [ opt.scVec( slew_idxs(1) ); opt.scVec( slew_idxs ) ];

slew_i( 1:opt.varAmts( opt.varNumDictionary( varName ) ) ) =...
    ( 1:opt.varAmts( opt.varNumDictionary( varName ) ) ).';
slew_j( 1:opt.varAmts( opt.varNumDictionary( varName ) ) ) =...
    slew_idxs;
slew_v( 1:opt.varAmts( opt.varNumDictionary( varName ) ) ) =...
    ones( size( slew_idxs ) );

% -1 on first subdiagonal
slew_i( ( opt.varAmts( opt.varNumDictionary( varName ) ) + 1 ):end ) =...
    ( 2: ( opt.varAmts( opt.varNumDictionary( varName ) ) + 1 ) ).';
slew_j( ( opt.varAmts( opt.varNumDictionary( varName ) ) + 1 ):end ) =...
    slew_idxs;
slew_v( ( opt.varAmts( opt.varNumDictionary( varName ) ) + 1 ):end ) =...
    -1 * ones( size( slew_idxs ) );

% Generate matrices
A_slew_pos = sparse( slew_i, slew_j, slew_v,...
    opt.varAmts( opt.varNumDictionary( varName ) ) + 1, opt.numVars );
% Generate "negative" difference matrix for breal_ORSP
A_slew_neg = -1 * A_slew_pos;

% Concatenate matrices
A_slew = [ A_slew_pos; A_slew_neg ];
b_slew = ( ( opt.dt * slewConstr ) ./ [ sc_vals; sc_vals] ) .*...
    ones( 2 * ( opt.varAmts( opt.varNumDictionary( varName ) ) + 1 ), 1 );
end