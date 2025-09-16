function [ NgramAssoc, NgramAmt, NgramDist ] = determineNgramRegions(...
    roi_body, rho_body, Ngrams, X, Y, Z, distnorm )
% This function will determine the Ngram association arrays. We will be
% using linear indexing for the size of the roi_body array.
%
% roi_body is logical array specifying which voxels are within the body
% rho_body is array of densities in kg/m^3
% Ngrams specifies how large the encompassing block of tissue should be (in
% terms of grams)
% X, Y, Z are the x, y, z positions of the voxels in ndgrid format
%
% Will return NgramAssoc matrix that will be a NumBodyVox x 1 cell array
% and each element of the cell array will have a vector of linear indices
% telling the user which voxels are within the Ng block of tissue
arguments
    roi_body
    rho_body
    Ngrams
    X
    Y
    Z
    distnorm = 2 % the norm to calculate distance from the target voxel
end

%% Determine domain size 
domSize = size( roi_body );

%% Determine linear indices for the roi_body indices
[ I, J, K ] = ndgrid(...
    1:domSize( 1 ),...
    1:domSize( 2 ),...
    1:domSize( 3 )...
    );

ind_3D = zeros( domSize );
ind_3D( : ) = 1:prod( domSize );

%% calculate dV
x = squeeze( X( :, 1, 1 ) );
y = squeeze( Y( 1, :, 1 ) );
z = squeeze( Z( 1, 1, : ) );

dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);

dV = dx * dy * dz;

%% Calculate universal offset points
offsetAmt = 0.036;
xo = -offsetAmt : dx : offsetAmt;
yo = -offsetAmt : dy : offsetAmt;
zo = -offsetAmt : dz : offsetAmt;

[ Xo, Yo, Zo ] = ndgrid( xo, yo, zo );

ro = [ Xo( : ), Yo( : ), Zo( : ) ];
norm_ro = vecnorm( ro, distnorm, 2 );

[ Io, Jo, Ko ] = ndgrid( 1:length(xo), 1:length(yo), 1:length(zo) );

% Sort based on magnitude of distance from point of interest
[ norm_ro, sort_norm_ro ] = sort( norm_ro, 1, 'ascend' );

epsround = eps( 1e1 );
ijkoi = [ find( abs( xo ) < epsround ), find( abs( yo ) < epsround ), find( abs( zo ) < epsround) ]; % find points that are at ro = [ 0, 0, 0 ]
ijko = [ Io( : ), Jo( : ), Ko( : ) ];
ijko =  ijkoi - ijko( sort_norm_ro, : );

%% Create vector of spatial positions and rho
ijk = [ I( roi_body ), J( roi_body ), K( roi_body ) ];
numPointsFOV = size( ijk, 1 );
mass_body = dV * rho_body;

%% Determine Nkg constraint
Nkg = Ngrams / 1e3; % convert grams to kg

%% Iterate over all points in the roi_body

% Initialize output arrays
NgramAssoc = cell( numPointsFOV, 1 );
NgramAmt = zeros( numPointsFOV, 1 );
NgramDist = zeros( numPointsFOV, 1 );

for ii = 1:numPointsFOV % iterate over each point
% for ii = 1:1000 % iterate over each point
    
    % Determine ijk for this point
    ijk_i = ijk( ii, : );
    
    % Determine ball of shortest norm around point of interest
    ijk_sort_normr = ijk_i + ijko;

    % Ensure that indices are within domain bounds
    ijk_sort_lb = all( ijk_sort_normr >= 1, 2 );
    ijk_sort_ub = all(...
        [...
        ( ijk_sort_normr(:, 1) <= domSize( 1 ) ),...
        ( ijk_sort_normr(:, 2) <= domSize( 2 ) ),...
        ( ijk_sort_normr(:, 3) <= domSize( 3 ) ) ], 2 );
    ijk_sort_valid = ijk_sort_lb & ijk_sort_ub;
    ijk_sort_normr = ijk_sort_normr( ijk_sort_valid, : );

    % lini_sort_normr = sub2ind(...
    %     domSize, ijk_sort_normr( :, 1 ), ijk_sort_normr( :, 2 ), ijk_sort_normr( :, 3 ) );
    % Following code is significantly faster than sub2ind
    lini_sort_normr = ...
        ijk_sort_normr( :, 1 ) +...
        ( ijk_sort_normr( :, 2 ) - 1 ) * domSize( 1 ) +...
        ( ijk_sort_normr( :, 3 ) - 1 ) * domSize( 1 ) * domSize( 2 );
    
    % Only use points that are within the body
    roi_body_sort_normr = roi_body( lini_sort_normr );
    lini_sort_normr = lini_sort_normr( roi_body_sort_normr );
    
    % determine sorted mass and indices around point of interest
    mass_sort_normr = mass_body( lini_sort_normr );
    ind_sort_normr = ind_3D( lini_sort_normr );
    
    % Determine cumulative mass moving away from point of interest in terms of
    % whatever norm we selected above
    massicumsum = cumsum( mass_sort_normr );
    
    % Determine first point that meets threshold
    cumsumidxthresh = find( massicumsum >= Nkg, 1, 'first' );

    NgramAssoc{ ii } = uint32( ind_sort_normr( 1:cumsumidxthresh ) );
    NgramAmt( ii ) = massicumsum( cumsumidxthresh );
    NgramDist( ii ) = norm_ro( cumsumidxthresh );
    
end

end