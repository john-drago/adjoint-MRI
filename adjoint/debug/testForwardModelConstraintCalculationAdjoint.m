function [ c, gradc, adjointToc ] = testForwardModelConstraintCalculationAdjoint( pSc, opt, binProfile, numLoops )

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
    [ c, ~, gradc, ~ ] = opt.nonlcon( pSc );
end
adjointToc = toc( adjointTic ) ;
fprintf('Analytic Constraint Gradient Time:\t %.5g\n', adjointToc);
if binProfile
    profile off;
    profile viewer;
end

end