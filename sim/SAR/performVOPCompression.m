function [ VOPs, VOPtoc, VOPID ] = performVOPCompression( QLocalMats, algorithm, varargin )
% This function will take in a set of QLocal matrices that determine the
% SAR in many locations from a given RF input. For example:
%
% SAR(r) = x^H Q x
%
% With this set of QLocal matrices, we can then determine VOP matrices that
% place an upperbound on the peak localSAR for all of the considered points
% in the domain
%
% QLocalMats: numRFChannels x numRFChannels x numFOVVoxels
% algorithm: specifies the VOP compression algorithm, options:
%   - eichfelder-fixed-overestimation
%   - eichfelder-determine-overestimation
%   - lee
%   - kuehne
%

%% Run VOP compression
VOPtic = tic;
if length( size( QLocalMats ) ) == 3

    if matches( algorithm, "eichfelder-fixed-overestimation", "IgnoreCase", true )
        error( "Not implemented yet" )
    elseif matches( algorithm, "eichfelder-determine-overestimation", "IgnoreCase", true )
        error( "Not implemented yet" )
    elseif matches( algorithm, "lee", "IgnoreCase", true )
        % Lee, MRM, 2012. doi: 10.1002/mrm.23140
        %
        % For the Lee algorithm, we assume that varargin contains the following
        % parameters:
        %   - epsG
        %   - overestMat
        %
        numRFchannels = size( QLocalMats, 1 );
        [ epsG, ovMat ] = processVarargin_LeeKuehne( varargin, numRFchannels );

        [ VOPs, VOPID ] = performVOPCompression_Lee( QLocalMats, epsG, ovMat );

    elseif matches( algorithm, "kuehne", "IgnoreCase", true )
        % Kuehne, ISMRM #478, 2017.
        %
        % For the Kuehne algorithm, we assume that varargin contains the following
        % parameters:
        %   - epsG
        %   - overestMat
        %
        numRFchannels = size( QLocalMats, 1 );
        [ epsG, ovMat ] = processVarargin_LeeKuehne( varargin, numRFchannels );

        [ VOPs, VOPID ] = performVOPCompression_Kuehne( QLocalMats, epsG, ovMat );
    else
        error( "Unknown compression algorithm" )
    end

else
    [ VOPs, VOPID ] = max( QLocalMats );
end
VOPtoc = toc( VOPtic );

end

%% Helper Function
% ----------------------------------------------------------------------- %
function [ epsG, ovMat ] = processVarargin_LeeKuehne( vararginCell, numRFchannels )

vararginlength = length( vararginCell );

% Overestimation factor
if vararginlength > 0
    epsG = vararginCell{ 1 }; % set the overestimation factor
else
    epsG = 0.10; % default 10% overestimation
end

% Set the overestimation matrix
if vararginlength > 1
    overestMat = vararginCell{ 2 }; % set the overestimation factor

    % if user supplied an overestimation matrix, then use that.
    % Otherwise, we can determine what to do from use input
    if isnumeric( overestMat )
        overestMatSz = size( overestMat );
        if ( overestMatSz(1) == size( QLocalMats, 1 ) ) & ( overestMatSz(2) == size( QLocalMats, 2 )  )
            ovMat = overestMat;
        elseif isscalar( overestMat )
            ovMat = overestMat * eye( size( QLocalMats, 1 ) );
        end

    elseif ischar( overestMat ) | isstring( overestMat )
        overestMat = char( overestMat );
        if matches( overestMat, 'identity', 'IgnoreCase', true )
            ovMat = eye( size( QLocalMats, 1 ) );
        elseif matches( overestMat, 'maxlocalSAReig', 'IgnoreCase', true )
            ovMat = 'maxlocalSAReig';
        end
    end

else
    ovMat = eye( numRFchannels );
end

end
% ----------------------------------------------------------------------- %