function [ brainMask, spatialWeight ] = generateMaskFromB0(...
    db0_mns, body_mask, b1p_CP_mns )

gyro = 267.5e6; % define gyromagnetic ratio

% define dummy dTE to generate phase distribution to then calculate weight
% metrics from the ROMEO paper with
dummydTE = 3.06e-3; 

% phony phase difference
phasediff = wrapToPi( db0_mns * gyro * dummydTE );

% calculate spatial weights;
spatialWeightOrig = calcSpatialWeight( phasediff, body_mask );

spatialWeight = spatialWeightOrig .* abs( b1p_CP_mns );
spatialWeight = spatialWeight / max( spatialWeight, [], 'all' );

% load example nifti image
brainMask = generateBrainMaskBET( spatialWeight );

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function brainMask = generateBrainMaskBET( spatialWeight )
currDir = fileparts( mfilename( 'fullpath' ) );
% currDir_linux = replace( currDir, ' ', '\ ' );

% specify dummy NIFTI file
mask_path = fullfile( currDir, 'dummy_nifti.nii' );
mask_path_linux = replace( mask_path, ' ', '\ ' );
mask_nii = load_nifti( mask_path );

mask_path_bet = fullfile( strcat( mask_path(1:end-4), "_BET.nii" ) );
mask_path_bet = char(mask_path_bet);
mask_path_bet_linux = replace( mask_path_bet, ' ', '\ ' );

mask_path_bet_gz = fullfile( strcat( mask_path_bet_linux, ".gz" ) );
mask_path_bet_gz = char(mask_path_bet_gz);
mask_path_bet_gz_linux = replace( mask_path_bet_gz, ' ', '\ ' );

if isfile( mask_path_bet_linux )
    [ ~, ~ ] = system( sprintf( "rm %s", mask_path_bet_linux ) );
end

if isfile( mask_path_bet_gz_linux )
    [ ~, ~ ] = system( sprintf( "rm %s", mask_path_bet_gz_linux ) );
end

% reconfigure nifit file
mask_nii.vol = spatialWeight; % add the spatial weight mask
mask_nii.dim( 1 ) = 3; % first entry is the number of dimensions
mask_nii.dim( 2:4 ) = size( spatialWeight ); % set spatial dimensions
mask_nii.dim( 5:end ) = 1; % make everything else 1
mask_nii.pixdim( 2 ) = 2;
mask_nii.pixdim( 3 ) = 2;
mask_nii.pixdim( 4 ) = 2;
mask_nii.pixdim( 5 ) = 1; % set time scale to be 1
mask_nii.xyzt_units = 2; % using millimeters binary
mask_nii.datatype = 64; % code for double entries
mask_nii.bitpix = 64; % number of bits per pixel
% mask_nii.quatern_x = 0; % unknown transformation
% mask_nii.quatern_y = 0; % unknown transformation
% mask_nii.quatern_z = 0; % unknown transformation
mask_nii.qform_code = 0; % unknown transformation
mask_nii.sform_code = 0; % unknown transformation

% save nifti file
save_nifti( mask_nii, mask_path );

% use BET to generate brain mask
procPath = fileparts( mask_path );

compType = computer;
if ~strcmpi( compType(1:4), "glnx" )
    error( "Unknown how to process in non-Linux system." )
end

% Generate temp file
tempBETFile = fullfile( currDir, "temp_bet.sh" );
tempBETFile_linux = replace( tempBETFile, ' ', '\ ' );
tempFileID = fopen( tempBETFile, 'w' );

fprintf( tempFileID, "#!/bin/bash\n" );
fprintf( tempFileID, "\n" );
fprintf( tempFileID, "source /usr/local/freesurfer/nmr-dev-env-bash\n" );
fprintf( tempFileID, "\n" );
fprintf( tempFileID, "bet %s %s -f 0.8125 -g +0.000 -R \n", mask_path_linux, mask_path_bet_linux );
fprintf( tempFileID, "\n" );
fprintf( tempFileID, "gunzip %s\n", strcat(mask_path_bet_linux, ".gz" ) );
fprintf( tempFileID, "\n" );

fclose( tempFileID );

[ ~, ~ ] = system( sprintf( "chmod 750 %s", tempBETFile_linux) );
% [ ~, cmdoutChmodBET ] = system( sprintf( "chmod 750 %s", tempBETFile_linux) );

[ status, ~ ] = system( tempBETFile );
% [ status, cmdoutTempBET ] = system( tempBETFile_linux );

if status == 0
    fprintf( "\n" );
    fprintf( "Brain extraction tool (BET) successful in command line.\n" );
    fprintf( "%s\n", procPath );
    fprintf( "\n" );
else
    error( "Brain extraction tool (BET) unsuccessful.\n" );
end

[ ~, ~ ] = system( sprintf( "rm %s", tempBETFile_linux ) );
% [ status, cmdoutDeleteBETTemp ] = system( sprintf( "rm %s", tempBETFile ) );

% Now load 

brain_mask_nii = load_nifti( mask_path_bet );
brainMask = ( brain_mask_nii.vol > 0 );

if isfile( mask_path_bet_gz_linux )
    [ ~, ~ ] = system( sprintf( "rm %s", mask_path_bet_gz_linux ) );
end

if isfile( mask_path_bet_linux )
    [ ~, ~ ] = system( sprintf( "rm %s", mask_path_bet_linux ) );
end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function spatialWeight = calcSpatialWeight( phasediff, body_mask )

sz_phasediff = size( phasediff );
[ Io, Jo, Ko ] = ndgrid(...
    1 : sz_phasediff(1),...
    1 : sz_phasediff(2),...
    1 : sz_phasediff(3) );
ijko = [ Io( : ), Jo( : ), Ko( : ) ];

expand_phasediff = zeros( ( size( phasediff ) + 2 ) );
expand_phasediff( (2:(end-1)), (2:(end-1)), (2:(end-1)) ) = phasediff;

sz_expand_phasediff = size( expand_phasediff );

[ Ie, Je, Ke ] = ndgrid(...
    1 : sz_expand_phasediff(1),...
    1 : sz_expand_phasediff(2),...
    1 : sz_expand_phasediff(3) );

Ie_base = Ie(...
    (2:(sz_expand_phasediff(1)-1)),...
    (2:(sz_expand_phasediff(2)-1)),...
    (2:(sz_expand_phasediff(3)-1)) );
Je_base = Je(...
    (2:(sz_expand_phasediff(1)-1)),...
    (2:(sz_expand_phasediff(2)-1)),...
    (2:(sz_expand_phasediff(3)-1)) );
Ke_base = Ke(...
    (2:(sz_expand_phasediff(1)-1)),...
    (2:(sz_expand_phasediff(2)-1)),...
    (2:(sz_expand_phasediff(3)-1)) );
ijke_base = [ Ie_base( : ), Je_base( : ), Ke_base( : ) ];
lind_base = sub2ind( sz_expand_phasediff, ijke_base(:,1), ijke_base(:,2), ijke_base(:,3) ); 

spatialWeightArr = zeros( [ size(phasediff), 6 ] );

swidx = 0;
for ax = 1:3
    for dr = [ -1, 1 ]

        shift_sub = zeros( 1, 3 );
        shift_sub( ax ) = dr;
        ijke_shift = ijke_base + shift_sub;
        lind_shift = sub2ind( sz_expand_phasediff, ijke_shift(:,1), ijke_shift(:,2), ijke_shift(:,3) ); 
        
        weight_ax_dr = 1 - abs( wrapToPi( expand_phasediff( lind_shift ) - expand_phasediff( lind_base ) ) / pi );
        
        swidx = swidx + 1;
        
        lind_sw = sub2ind( [ sz_phasediff, 6 ], ijko(:,1), ijko(:,2), ijko(:,3), swidx*ones(size(ijko(:,1))) );

        spatialWeightArr( lind_sw ) = weight_ax_dr;
    end
end

spatialWeight = prod( spatialWeightArr, 4 );
spatialWeight = spatialWeight .* body_mask;
end
% ----------------------------------------------------------------------- %