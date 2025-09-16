function [ QMat ] = determineGlobalSARMatrix( em )
% This function computes the Q matrix average over all of the indices in
% the body using the fields in the em struct.
%
% em is a structure that contains: X, Y, Z, roi_body, E, rho_body,
% sigma_body for SAR calculation
%
% We will output QMat which will be a numChannel x numChannel matrix
% array that contains information for calculating global SAR. 
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

%% Get the linear indices for all of the voxels in the field of view
domSize = size( roi_body );
ind_3D = zeros( domSize );
ind_3D( : ) = 1:prod( domSize );
linind = ind_3D( roi_body );

%% Reshape E field for faster indexing
% domSize = size( roi_body );
E = reshape( E, numRF, numel( roi_body ), 3 ); % reshape to the size of the domain

% calculate coefficient to E^H E before the loop (some values will be NaNs)
rhoVec = rho_body( : );
sigmaVec = sigma_body( : );
coeffVec = 0.5 * ( sigmaVec ./ rhoVec );
% still will need to normalize by ΔV / V (which is 1 / numSumVox, assuming
% the voxels are all the same size)

%% Iterate over each of the voxels and sum

numSumVox = length( linind ); % determine number of voxels in summation

% get E field for the voxels in summation in form:
% 3 x Num RF Channels x num Vox
E_voxels = permute( E( :, linind, : ), [ 3 1 2 ] );

% perform the E^H E operation pagewise and then scale with the coeff vector
Qvox = reshape( coeffVec( linind ), [ 1 1 numSumVox ] ) .* ...
    pagemtimes( E_voxels, 'ctranspose', E_voxels, 'none' );

% sum over the voxels for this Ngram group and scale by voxel size
QMat = 1/numSumVox * sum( Qvox, 3 );

end