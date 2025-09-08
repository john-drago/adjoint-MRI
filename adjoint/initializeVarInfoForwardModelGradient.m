function varInfo = initializeVarInfoForwardModelGradient( opt )

varInfo = struct;
varInfo.varCtr = uint32( 0 );
varInfo.varAmtsCtr = uint32( 0 );
varInfo.varAmts = zeros( opt.estMaxActiveVarsTimeStep, 1, "uint32" );
varInfo.varIdxsStartStop = zeros( opt.estMaxActiveVarsTimeStep, 2, "uint32" );
varInfo.varIdxs = zeros( opt.estMaxActiveVarsTimeStep, 1, "uint32" );

end