function [ A_slew, b_slew ] = constraintvarkTPPulseSlew(...
    opt, varname, slewConst )
% Function will ensure that the slew constraints are honored for the
% variables considered

% Get number kT-points
num_kTP = opt.num_kTP;

% get variable positions
var_idx = opt.( strcat(varname, '_idx') );
dt_idx = opt.dt_idx;
blip_dt_idx = dt_idx( 2:2:end );

% determine number of variables
numSlewVarPerkTP = size( var_idx, 1 ) / ( num_kTP-1 );
var_idx_rshp = reshape( var_idx, [ numSlewVarPerkTP, (num_kTP-1) ] );

% initialize matrix
A_slew = zeros( 2 * numSlewVarPerkTP * ( num_kTP-1 ) , opt.numVars );

for nn = 1:( num_kTP-1 )
    A_row_pos = ( ( 2 * numSlewVarPerkTP ) * ( nn - 1 ) + ( 1:numSlewVarPerkTP ) ).';
    A_row_neg = ( ( 2 * numSlewVarPerkTP ) * ( nn - 1 ) +...
        ( (numSlewVarPerkTP + 1):(2 * numSlewVarPerkTP) ) ).';

    pos_var_idx = sub2ind( size( A_slew ), A_row_pos, var_idx_rshp( :, nn ) );
    pos_dt_idx = sub2ind( size( A_slew ), A_row_pos, repmat( blip_dt_idx( nn ), [numSlewVarPerkTP 1] ) );
    neg_var_idx = sub2ind( size( A_slew ), A_row_neg, var_idx_rshp( :, nn ) );
    neg_dt_idx = sub2ind( size( A_slew ), A_row_neg, repmat( blip_dt_idx( nn ), [numSlewVarPerkTP 1] ) );
    
    % Positive slew
    % magnitude of blip
    A_slew( pos_var_idx ) = opt.scVec( var_idx_rshp( :, nn ) ).';
    % dt length
    A_slew( pos_dt_idx ) =...
        repmat( - opt.scVec( blip_dt_idx( nn ) ) * slewConst / 2, [ numSlewVarPerkTP, 1 ] );

    % Negative slew
    % magnitude of blip
    A_slew( neg_var_idx ) = - opt.scVec( var_idx_rshp( :, nn ) ).';
    % dt length
    A_slew( neg_dt_idx ) =...
        repmat( - opt.scVec( blip_dt_idx( nn ) ) * slewConst / 2, [ numSlewVarPerkTP, 1 ] );
end
b_slew = zeros( 2 * numSlewVarPerkTP * ( opt.num_kTP-1 ), 1 );

end