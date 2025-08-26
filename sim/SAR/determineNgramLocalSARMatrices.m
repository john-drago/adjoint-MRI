function [ QMats, QMatsIdx ] = determineNgramLocalSARMatrices( NgramAssoc, em )
% This function takes in a set of Ngram association arrays and computes the
% Q matrix average over all of the indices that constitute an associate
% array according to the fields in the em struct.
%
% NgramAssoc is the array that contains the linear indices for the 3D array
% that specify the Ngram regions aroung the first linear index.
%
% em is a structure that contains: X, Y, Z, roi_body, E, rho_body,
% sigma_body for SAR calculation
%
% We will output QMats which will be a numChannel x numChannel x numVox
% array that contains information for calculating Ng Local SAR. 
% QMatInd specifies the linear index for the location of the Ng Local SAR Q
% matrix
% 

%% Extract from the em struct
roi_body = em.roi_body;
rho_body = em.rho_body;
sigma_body = em.sigma_body;
E = em.E;
% X = em.X;
% Y = em.Y;
% Z = em.Z;

%% Determine Number of Coils
Esize = size( E );
if ( length( Esize ) == 4 ) & ( Esize( end ) == 3 ) % if the E array contains the fields for a single-channel transmit coil
    numRF = 1;
elseif ( length( Esize ) == 5 ) & ( Esize( end ) == 3 ) % if the E array contains the fields for a multichannel transmit coil
    numRF = Esize( 1 );
end

%% calculate dV
% x = squeeze( X( :, 1, 1 ) );
% y = squeeze( Y( 1, :, 1 ) );
% z = squeeze( Z( 1, 1, : ) );
% 
% dx = x(2) - x(1);
% dy = y(2) - y(1);
% dz = z(2) - z(1);

% dV = dx * dy * dz;

%% Initialize Arrays
QMats = zeros( numRF, numRF, length( NgramAssoc ) );
QMatsIdx = zeros( length( NgramAssoc ), 1 );

%% Reshape E field for faster indexing
% domSize = size( roi_body );
E = reshape( E, numRF, numel( roi_body ), 3 ); % reshape to the size of the domain

% calculate coefficient to E^H E before the loop (some values will be NaNs)
rhoVec = rho_body( : );
sigmaVec = sigma_body( : );
coeffVec = 0.5 * ( sigmaVec ./ rhoVec );
% still will need to normalize by ΔV / V (which is 1 / numSumVox, assuming
% the voxels are all the same size)

%% Iterate over each of the voxels
for nn = 1 : length( NgramAssoc )
    
    assoc_i = NgramAssoc{ nn }; % get linear indices for Ngram average
    numSumVox = length( assoc_i ); % determine number of voxels in summation
    
    % get E field for the voxels in summation in form: 
    % 3 x Num RF Channels x num Vox
    E_voxels = permute( E( :, assoc_i, : ), [ 3 1 2 ] );
    
    % perform the E^H E operation pagewise and then scale with the coeff vector
    Qvox = reshape( coeffVec( assoc_i ), [ 1 1 numSumVox ] ) .* ...
        pagemtimes( E_voxels, 'ctranspose', E_voxels, 'none' );
    
    % sum over the voxels for this Ngram group and scale by voxel size
    QMats( :, :, nn ) = 1/numSumVox * sum( Qvox, 3 );
    
    % get index for where this voxel is
    QMatsIdx( nn ) = assoc_i( 1 );

end

%% Post-process
% if there is a single-channel then make Q matrix into a column vector
if numRF == 1
    QMats = squeeze( QMats );
end

end