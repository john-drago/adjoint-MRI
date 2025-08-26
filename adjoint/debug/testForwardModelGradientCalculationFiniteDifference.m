function [ J_FD, gradp_J_FD, fdToc ] = testForwardModelGradientCalculationFiniteDifference( pSc, opt )

%% Test analytic
costfn = @( pSc ) opt.runCostFunction( pSc, opt );

fdTic = tic;
[ dJdp_FD, J_FD ] = calcJacobianFiniteDifference( costfn, pSc, 1e-7 );
fdToc = toc( fdTic );
fprintf('FD Cost Gradient Time:\t\t %.5g sec\n', fdToc);
gradp_J_FD = dJdp_FD.';

end