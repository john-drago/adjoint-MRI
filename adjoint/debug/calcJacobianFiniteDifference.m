function [ dfdx, f_x0 ] = calcJacobianFiniteDifference( evalf, x0, dxFD )
% This function will calculate the finite difference jacobian
arguments
    evalf; % function to compute Jacobian with
    x0; % point at which to approximate Jacobian
    dxFD = 1e-6; % finite difference step size
end

x0 = x0( : );
xk = x0;
N = length( x0 );

f_x0 = evalf( x0 );
dfdx = zeros( length( f_x0 ), N );

for nn = 1:N
    xk( nn ) = xk( nn ) + dxFD;
    f_xk = evalf( xk );
    dfdx( :, nn ) = ( f_xk - f_x0 ) / dxFD;
    xk( nn ) = x0( nn );
end

end