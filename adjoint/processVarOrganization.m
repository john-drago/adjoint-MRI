function [ opt, varArray, varNames, varAmts, varNumDictionary, varCumSum, varIdxs, numVars ] =...
    processVarOrganization( opt, varArray )

varNumDictionary = dictionary();

varCumSum = cumsum( [ 1; cell2mat(varArray( :, 2 ) ) ] );
varIdxs = cell( (size( varCumSum, 1 ) - 1), 1 );
for ii = 1 : size( varIdxs, 1 )
    varIdxs{ ii } = uint32( varCumSum( ii ) : ( varCumSum( ii+1 ) - 1 ) ).';
    varNumDictionary( varArray{ ii, 1 } ) = ii;
end
varCumSum = uint32( varCumSum( 1 : (end-1) ) );

varArray = [ varArray, varIdxs ];
varNames = string( varArray( :, 1 ) );
varAmts = uint32( cell2mat( varArray( :, 2 ) ) );

numVars = varCumSum( end ) + varAmts( end ) - 1;

varCumSum = cumsum( varAmts );

opt.varNames = varNames;
opt.varAmts = varAmts;
opt.varNumDictionary = varNumDictionary;
opt.varCumSum = varCumSum;
opt.varIdxs = varIdxs;
opt.numVars = numVars;

for vv = 1:length( varNames )
    opt.( sprintf( "%s_idx", replace(varNames( vv ), "-", "_" ) ) ) = varIdxs{ vv };
end

end