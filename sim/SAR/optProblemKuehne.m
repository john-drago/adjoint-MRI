function [ c_wv_opt, out ] = optProblemKuehne( c_wv_0, Qcurrent, Qsub, ovMat, algorithm )
arguments
    c_wv_0
    Qcurrent
    Qsub
    ovMat
    algorithm = 1; % 1: uses Orzada implementation........ 2: use Kuehne original implementation
end

fminopt = optimoptions( 'fmincon' );

if algorithm == 1
    fminopt.Algorithm = 'sqp';
elseif algorithm == 2
    fminopt.Algorithm = 'interior-point';
else
    error( "Unknown Kuehne algorithm implementation" )
end

tol = 1e-14;

fminopt.EnableFeasibilityMode = true;
fminopt.ConstraintTolerance = tol;
fminopt.MaxFunctionEvaluations = inf;
fminopt.OptimalityTolerance = tol;
fminopt.StepTolerance = tol;

fminopt.Display = 'off'; % 'notify-detailed';
fminopt.MaxIterations = 4e2;

Aeq = [];
beq = [];

lb = zeros( size( c_wv_0 ) );
ub = ones( size( c_wv_0 ) );

A = ones( 1, size( c_wv_0, 1 ) );
b = 1;

if algorithm == 1
    fminopt.SpecifyObjectiveGradient = false;
    objfun = @( c_wv ) objCostKuehne1( c_wv(:), Qcurrent, Qsub, ovMat );
    nonlcon = [];
elseif algorithm == 2
    fminopt.SpecifyObjectiveGradient = true;
    grad = ones( size( c_wv_0(:) ) );
    objfun = @( c_wv ) objCostKuehne2( c_wv(:), grad );
    nonlcon = @( c_wv ) PSDConstraint( c_wv, Qcurrent, Qsub, ovMat );
end

[ c_wv_opt, fval, exitflag, output ] = fmincon(...
    objfun,...
    c_wv_0, A, b, Aeq, beq, lb, ub, nonlcon, fminopt );

if algorithm == 1
    out = fval;
elseif algorithm == 2
    if exitflag < 1
        out = 1;
    else
        out = -1;
    end
    
end

end
%% Helper Function
% ----------------------------------------------------------------------- %
function [ output ] = objCostKuehne1( c_wv, Qcurrent, Qsub, ovMat )
% Orzada implementation
c_wv = reshape( c_wv, [ 1, 1, length(c_wv) ] );
% LMI = sum( ( Qsub ) .* c_wv, 3 ) + ovMat - Qcurrent;
LMI = sum( ( Qsub + ovMat ) .* c_wv, 3 ) - Qcurrent;
eigLMI = eig( LMI );
mineig = min( real( eigLMI ) );
output = - mineig;
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ output, grad ] = objCostKuehne2( c_wv_opt, grad )
% Kuehne original implementation
output = sum( c_wv_opt );
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ c, ceq ] = PSDConstraint( c_wv, Qcurrent, Qsub, ovMat )
% Kuehne original implementation
ceq = [];

c_wv = reshape( c_wv, [ 1, 1, length(c_wv) ] );
LMI = sum( ( Qsub + ovMat ) .* c_wv, 3 ) - Qcurrent;
eigLMI = eig( LMI );
mineig = min( real( eigLMI ) );
c = - mineig;

end
% ----------------------------------------------------------------------- %