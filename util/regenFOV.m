function Array3D = regenFOV(vec, roiFOV)
% Function that will regenerate a region of interest within the FOV if it
% is given a structure that has the roiFOV (logical) array as a field.
arguments
    vec % vector to put into a 3D array
    roiFOV % pass in specific roiFOV that you want to be used
end

% Make sure vector is a column vec
vec = vec(:);

% Check to make sure that roiFOV and vector are compatible 
if sum(roiFOV, 'all') ~= size(vec,1)
    error("Incorrect roiFOV indexing for vector")
end

% Initialize Array
Array3D = zeros(size(roiFOV));

% Create indexing based on the roiFOV logical entries
[II, JJ, KK] = ndgrid(1:size(roiFOV,1), 1:size(roiFOV,2), 1:size(roiFOV,3));
ijk = [...
    II(logical(roiFOV)),...
    JJ(logical(roiFOV)),...
    KK(logical(roiFOV))];
if length(size(ijk))==3
    ijk = squeeze(ijk)';
end

% Put into 3D array
indArray = sub2ind(size(roiFOV), ijk(:,1), ijk(:,2), ijk(:,3));
Array3D(indArray) = vec;
end