function [ c_wv_opt, converged ] = feasibilityProblemLee( c_wv_0, Qcurrent, Qsub, ovMat )

tol = 1e-14;

fminopt = optimoptions( 'fmincon' );
fminopt.Algorithm = 'interior-point';
fminopt.EnableFeasibilityMode = true;
fminopt.ConstraintTolerance = tol;
fminopt.MaxFunctionEvaluations = inf;
fminopt.OptimalityTolerance = tol;
% fminopt.Display = 'notify-detailed';
fminopt.Display = 'off';
fminopt.MaxIterations = 4e2;

Aeq = [];
beq = [];

lb = zeros( size( c_wv_0 ) );
ub = ones( size( c_wv_0 ) );

A = ones( 1, size( c_wv_0, 1 ) );
b = 1;

nonlcon = @( c_wv ) PSDConstraint( c_wv, Qcurrent, Qsub, ovMat );

[ c_wv_opt, ~, exitflag, ~ ] = fmincon(...
    @( c_wv_0 ) objFeasibility( c_wv_0(:) ),...
    c_wv_0, A, b, Aeq, beq, lb, ub, nonlcon, fminopt );

if exitflag ~= 1
    converged = false;
else
    converged = true;
end

end
%% Helper Function
% ----------------------------------------------------------------------- %
function output = objFeasibility( ~ )
output = 0;
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ c, ceq ] = PSDConstraint( c_wv, Qcurrent, Qsub, ovMat )
ceq = [];

c_wv = reshape( c_wv, [ 1, 1, length(c_wv) ] );
LMI = sum( ( Qsub ) .* c_wv, 3 ) + ovMat - Qcurrent;
eigLMI = eig( LMI );
mineig = min( real( eigLMI ) );
c = - mineig;

end
% ----------------------------------------------------------------------- %