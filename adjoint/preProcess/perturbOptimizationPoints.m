function [ x_pert, y_pert, z_pert ] =  perturbOptimizationPoints(...
    x_opt, y_opt, z_opt, x_const, y_const, z_const )

x_pert = perturbOptPoints( x_opt, x_const );
y_pert = perturbOptPoints( y_opt, y_const );
z_pert = perturbOptPoints( z_opt, z_const );

end

%% Helper Function
% ----------------------------------------------------------------------- %
function w_pert = perturbOptPoints( w_opt, w_const )

idxtol = 1e-8;

diffw1 = diff( w_opt( 1:2 ) );
diffwend = diff( w_opt( (end-1):end ) );
w_opt_extend = [ w_opt( 1 ) - diffw1, w_opt, w_opt( end ) + diffwend ];
w_opt_samp = mean(...
    [...
    w_opt_extend( 1:(end-1) );...
    w_opt_extend( 2:(end) ) ], 1 );
dw_opt = diff( w_opt_samp );

w_pert = w_opt_samp( 1:(end-1) ) + dw_opt .* rand( size( dw_opt ) );

if ~isempty( w_const )
    w_const_rep = repmat( w_const( : ), [ 1, length( w_opt ) ] );
    w_const_idx = ( min( abs( w_opt - w_const_rep ), [], 1 ) < idxtol );
    w_pert( w_const_idx ) = w_opt( w_const_idx );
end

end
% ----------------------------------------------------------------------- %