function partialfpartialp = placepartialfpartialpCalculation( partialfpartialp, contribution, varInfo )

idxs = varInfo.varIdxsStartStop(varInfo.varCtr, 1):varInfo.varIdxsStartStop(varInfo.varCtr, 2);
partialfpartialp( :, :, idxs ) = contribution;

% if opt.useGPU
%     partialfpartialp( :, :, idxs ) = contribution;
% else
%     varAmts = varInfo.varAmts(varInfo.varCtr);
%     partialfpartialp( idxs ) =...
%         squeeze(...
%         mat2cell( contribution, opt.numPos, 3, ones( varAmts, 1 )  )...
%         );
% end

end