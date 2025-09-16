function [ dfdx, f_x0 ] = calcJacobianFiniteDifferenceArray( evalf, x0, dxFD )
% This function will calculate the finite difference jacobian
arguments
    evalf; % function to compute Jacobian with
    x0; % point at which to approximate Jacobian
    dxFD = 1e-7; % finite difference step size
end


x0_sz = size( x0 );
if x0_sz( end ) == 1
    x0_sznot1 = x0_sz( 1:(end-1) );
else
    x0_sznot1 = x0_sz;
end
x0_numel = prod( x0_sznot1 );
x0_sub = cell( size( x0_sz ) );

xk = x0;

f_x0 = evalf( x0 );
f_x0_sz = size( f_x0 );
f_x0_numel = prod( f_x0_sz );
if ( f_x0_sz( end ) == 1 ) && ( f_x0_sz( 1 ) == 1 )
    f_x0_sznot1 = [];
elseif f_x0_sz( end ) == 1
    f_x0_sznot1 = f_x0_sz( 1:(end-1) );
else
    f_x0_sznot1 = f_x0_sz;
end
fsubidx = cell( size( f_x0_sz ) );
[ fsubidx{:} ] = ind2sub( f_x0_sz, ( 1:f_x0_numel ).' );

if isempty( f_x0_sznot1 )
    f_x0_numel = 1;
else
    f_x0_numel = prod( f_x0_sznot1 );
end

dfdx = zeros( [ f_x0_sznot1, x0_sz ] );
dfdx_sz = size( dfdx );

for nn = 1:x0_numel

    [x0_sub{:}] = ind2sub( x0_sz, nn );
    x0_sub_expand = cellfun( @(c) repmat( c, [ f_x0_numel, 1 ] ), x0_sub, 'uniformoutput', false );
    
    if isempty(f_x0_sznot1)
        dfdx_sub_expand = x0_sub_expand;
    else
        dfdx_sub_expand = [fsubidx, x0_sub_expand];
    end

    xk( nn ) = x0( nn ) + dxFD;
    f_xk = evalf( xk );
    df = ( f_xk - f_x0 ) / dxFD;
    dfdxrshp = reshape( df, [ f_x0_numel, 1 ] );
    
    dfdx_ind = sub2ind( dfdx_sz, dfdx_sub_expand{:} );
    
    dfdx( dfdx_ind ) = dfdxrshp;
    xk( nn ) = x0( nn );
end

end