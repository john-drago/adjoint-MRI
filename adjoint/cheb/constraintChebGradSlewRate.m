function [ c, gradc ] = constraintChebGradSlewRate( pSc, opt, maxExtrema )

if nargin < 3
    maxExtrema = 1;
else
    maxExtrema = min( [ maxExtrema, opt.orderCheb_grad - 1 ] );
    maxExtrema = max( [ 1, maxExtrema ] );
end

gradSlewRateConstr = opt.gradSlewRate_constr;

% Assume pSc is scaled to be between -1 and 1
pSc = pSc( : );

tdom = opt.tdom;

grad_idx = opt.grad_idx;
grad_sc = opt.scVec( grad_idx );
grad = grad_sc .* pSc( grad_idx );
grad_idx_rshp = reshape( grad_idx, [ opt.numCheb_grad, 3 ]  );
grad_sc_rshp = reshape( grad_sc, [ opt.numCheb_grad, 3 ] );
grad_rshp = reshape( grad, [ opt.numCheb_grad, 3 ] );
dgrad_rshp = opt.D_grad * grad_rshp;
[ dgrad_minmax, dgrad_minmax_pos ] = chebMinMax( dgrad_rshp, tdom, maxExtrema );
if iscell( dgrad_minmax )
    dgrad_minmax = cell2mat( dgrad_minmax );
    dgrad_minmax_pos = cell2mat( dgrad_minmax_pos );
end
[ dgrad_absmax, dgrad_absmax_rowpos ] = maxk( abs( dgrad_minmax ), maxExtrema, 1 );

c_unsc = reshape( dgrad_absmax, [ maxExtrema*3, 1 ] ) - gradSlewRateConstr;
c = c_unsc / gradSlewRateConstr;

if nargout > 1

    gradc_unsc = zeros( opt.numVars, maxExtrema * 3 );

    dgrad_absmax_linpos = reshape( sub2ind( size(dgrad_minmax), dgrad_absmax_rowpos, repmat( 1:3, [ maxExtrema, 1 ] ) ), [ maxExtrema*3, 1 ] );
    dgrad_absmax_pos = reshape( dgrad_minmax_pos( dgrad_absmax_linpos ), [ maxExtrema*3, 1 ] );

    T_n =  transpose(...
        evalChebClenshaw( dgrad_absmax_pos, eye( opt.numCheb_grad - 1 ), tdom ) * opt.D_grad );
    % T_n_vec = reshape( T_n, [ ( opt.numCheb_grad ) * 3, 1 ] );

    % dgrad_idx_rshp = grad_idx_rshp( 1:end, : );
    % dgrad_sc_rshp = grad_sc_rshp( 1:end, : );

    dgradc_linpos = sub2ind(...
        size( gradc_unsc ),...
        reshape( permute( repmat( grad_idx_rshp, [ 1, 1, maxExtrema ] ), [ 1, 3, 2 ] ), [ opt.numCheb_grad, maxExtrema*3 ] ),...
        repmat( uint32( 1:(3*maxExtrema) ), [ opt.numCheb_grad, 1 ] ) );

    dgrad_minmax_sign_vec = ...
        repmat(...
        transpose( sign( dgrad_minmax( dgrad_absmax_linpos ) ) ), [ opt.numCheb_grad, 1 ] );

    gradc_unsc( dgradc_linpos ) =...
        ( dgrad_minmax_sign_vec .* T_n ) .* reshape( permute( repmat( grad_sc_rshp, [ 1, 1, maxExtrema ] ), [ 1, 3, 2 ] ), [ opt.numCheb_grad, maxExtrema*3 ] );

    gradc = gradc_unsc / gradSlewRateConstr;

end

end


% function [ c, gradc ] = constraintChebGradSlewRate( pSc, opt )
% 
% gradSlewRateConstr = opt.gradSlewRate_constr;
% 
% numCheb = opt.numCheb_grad;
% numCoils = 3;
% 
% numTimePoints = opt.numTimePoints;
% Tn = opt.Tn( :, 1:(numCheb-1) );
% dTn = Tn * opt.D_grad;
% 
% % Assume pSc is scaled to be between -1 and 1
% pSc = pSc( : );
% 
% grad_idx = opt.grad_idx;
% grad_sc = opt.scVec( grad_idx );
% grad = grad_sc .* pSc( grad_idx );
% grad_idx_rshp = reshape( grad_idx, [ numCheb, numCoils ]  );
% grad_sc_rshp = reshape( grad_sc, [ numCheb, numCoils ] );
% grad_rshp = reshape( grad, [ numCheb, numCoils ] );
% dgrad_time_rshp = dTn * grad_rshp;
% dgrad_abs_rshp = abs( dgrad_time_rshp );
% 
% c_unsc = reshape( dgrad_abs_rshp, [ numTimePoints*numCoils, 1 ] ) - gradSlewRateConstr;
% c = c_unsc / gradSlewRateConstr;
% 
% if nargout > 1
% 
%     gradc_unsc = zeros( opt.numVars, numTimePoints * numCoils );
% 
%     dgrad_abs_rshp_sign = ( reshape( sign( dgrad_time_rshp ), [ numTimePoints, 1, numCoils ] ) .* ( dTn ) )...
%         .* reshape( grad_sc_rshp, [ 1, numCheb, numCoils ] );
% 
%     grad_idx_gradc = repmat( reshape( grad_idx_rshp, [ 1, numCheb, numCoils ] ), [ numTimePoints, 1, 1 ] );
% 
%     c_idx_gradc = permute( repmat( reshape( uint32( 1:(numTimePoints*numCoils) ), [ numTimePoints, numCoils ] ), [ 1, 1, numCheb ] ), [ 1, 3, 2 ] );
% 
%     gradc_dgrad_idx = sub2ind( size( gradc_unsc ), grad_idx_gradc, c_idx_gradc );
% 
%     gradc_unsc( gradc_dgrad_idx ) = dgrad_abs_rshp_sign;
% 
%     gradc = gradc_unsc / gradSlewRateConstr;
% 
% end
% 
% end