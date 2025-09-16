function [ J_A, gradp_J_A, adjointToc ] = testForwardModelGradientCalculationAdjoint( pSc, opt, binProfile, numLoops )

if nargin < 3
    binProfile = false;
end
if nargin < 4
    numLoops = 1;
end

%% Test analytic
if binProfile
    profile on;
end
adjointTic = tic;
for ii = 1:numLoops
    [ J_A, gradp_J_A ] = opt.runCostFunction( pSc, opt );
end
adjointToc = toc( adjointTic ) ;
fprintf('Analytic Cost Gradient Time:\t %.5g sec\n', adjointToc);
if binProfile
    profile off;
    profile viewer;
end

end