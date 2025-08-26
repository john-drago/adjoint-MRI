function [ c_FD, gradc_FD, fdToc ] = testForwardModelConstraintCalculationFiniteDifference( pSc, opt )

fdTic = tic;
[ dcdp_FD, c_FD ] = calcJacobianFiniteDifference( opt.nonlcon, pSc, 1e-7 );
fdToc = toc( fdTic );
fprintf('FD Constraint Gradient Time:\t\t %.5g sec\n', fdToc);
gradc_FD = dcdp_FD.';

end