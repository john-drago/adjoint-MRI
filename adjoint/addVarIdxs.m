function varInfo = addVarIdxs( opt_idxs, varInfo )

varInfo.varCtr = varInfo.varCtr + 1;
varInfo.varAmts( varInfo.varCtr ) = numel( opt_idxs );

startStop = [ varInfo.varAmtsCtr+1, varInfo.varAmtsCtr + varInfo.varAmts( varInfo.varCtr ) ];
varInfo.varIdxsStartStop( varInfo.varCtr, : ) = startStop;

idxs = (varInfo.varAmtsCtr+1):(varInfo.varAmtsCtr+varInfo.varAmts( varInfo.varCtr ));
varInfo.varIdxs( idxs ) = opt_idxs;

varInfo.varAmtsCtr = varInfo.varAmtsCtr + varInfo.varAmts( varInfo.varCtr );

end