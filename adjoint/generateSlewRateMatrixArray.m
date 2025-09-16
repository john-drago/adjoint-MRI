function [ A_slew, b_slew ] = generateSlewRateMatrixArray( varName, numVars, opt, slewConstr )

varAmts = opt.varAmts( opt.varNumDictionary( varName ) );
varAmts_nn = varAmts / numVars;

slew_idxs = opt.varIdxs( opt.varNumDictionary( varName ) );
slew_idxs = slew_idxs{ 1 };

numRows = numVars * ( varAmts_nn + 1 );

% Create "positive" difference matrix
slew_i = zeros( varAmts*2, 1 );
slew_j = zeros( varAmts*2, 1 );
slew_v = zeros( varAmts*2, 1 );
sc_vals = zeros( numRows, 1 );

for nn = 1:numVars
    
    slew_i_nn = zeros( varAmts_nn*2, 1 );
    slew_j_nn = zeros( varAmts_nn*2, 1 );
    slew_v_nn = zeros( varAmts_nn*2, 1 );

    slew_idxs_nn = slew_idxs( ( (nn - 1)*varAmts_nn + (1 : varAmts_nn) ) );
    
    % 1 on the diagonal
    slew_i_nn( 1:varAmts_nn ) =...
        ( (nn - 1)*(varAmts_nn+1) + (1 : varAmts_nn) ).';
    slew_j_nn( 1:varAmts_nn ) =...
        slew_idxs_nn;
    slew_v_nn( 1:varAmts_nn ) =...
        ones( varAmts_nn, 1 );

    % -1 on first subdiagonal
    slew_i_nn( ( varAmts_nn + 1 ):end ) =...
        ( (nn - 1)*(varAmts_nn+1) + (2 : ( varAmts_nn + 1 ) ) ).';
    slew_j_nn( ( varAmts_nn + 1 ):end ) =...
        slew_idxs_nn;
    slew_v_nn( ( varAmts_nn + 1 ):end ) =...
        -1 * ones( varAmts_nn, 1 );

    sc_vals_nn = [ opt.scVec( slew_idxs_nn( 1 )) ; opt.scVec( slew_idxs_nn ) ];
    
    % Assign position
    slew_i( (nn - 1)*2*varAmts_nn + (1:2*varAmts_nn).' ) = slew_i_nn;
    slew_j( (nn - 1)*2*varAmts_nn + (1:2*varAmts_nn).' ) = slew_j_nn;
    slew_v( (nn - 1)*2*varAmts_nn + (1:2*varAmts_nn).' ) = slew_v_nn;

    sc_vals( (nn - 1)*( varAmts_nn + 1 ) + ( 1:( varAmts_nn+1 ) ).' ) =...
        sc_vals_nn;

end

% Generate matrices
A_slew_pos = sparse( slew_i, slew_j, slew_v,...
    numRows, opt.numVars );
% Generate "negative" difference matrix for breal_ORSP
A_slew_neg = -1 * A_slew_pos;

% Concatenate matrices
A_slew = [ A_slew_pos; A_slew_neg ];
b_slew = ( ( opt.dt * slewConstr ) ./ [ sc_vals; sc_vals] ) .*...
    ones( 2 * numRows, 1 );

end