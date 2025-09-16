function [ opt, val ] = addOptInfoVal( opt, val )
% This function will add the correct variable indices to the val struct
% from the opt struct.

val.varNames = opt.varNames;
val.varAmts = opt.varAmts;
val.varNumDictionary = opt.varNumDictionary;
val.varCumSum = opt.varCumSum;
val.varIdxs = opt.varIdxs;
val.numVars = opt.numVars;

end