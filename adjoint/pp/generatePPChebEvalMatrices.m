function [ varsToTimepoints, varsToChebByPeriods, varsToShapeByPeriods,...
    TnShape, shapeFnChebCoeffs ] =...
    generatePPChebEvalMatrices( numPPShape, numPPPeriods, PPVarIdxs, PPVarSFIdxs, PPSFIdxs, PPIdxs, PPStartStop, tvec )

numTimePoints = length( tvec );
numVarsPerChannel = PPVarIdxs( end, 2 );

% Mapping between the variables and shape functions
varsToShapeByPeriods = zeros( numPPShape * numPPPeriods, numVarsPerChannel );
shapeByPeriodsToVar = zeros( numVarsPerChannel, numPPShape * numPPPeriods );

for nn = 1:numPPPeriods

    % Create Mapping from the variables to shape functions
    SFIdxs = PPSFIdxs( nn, 1 ):PPSFIdxs( nn, 2 );
    numShapeNumPIdx = SFIdxs( PPVarSFIdxs( nn,1 ):PPVarSFIdxs( nn,2 ) );
    numVarsIdx = PPVarIdxs( nn,1 ):PPVarIdxs( nn,2 );

    varsToShapeByPeriodsLinIdx = ...
        sub2ind( size( varsToShapeByPeriods ), numShapeNumPIdx, numVarsIdx );

    varsToShapeByPeriods(...
        varsToShapeByPeriodsLinIdx ) = 1;

    % Create Mapping from the shape functions to variables
    shapeByPeriodsToVarsLinIdx = ...
        sub2ind( size( shapeByPeriodsToVar ), numVarsIdx, numShapeNumPIdx );

    shapeByPeriodsToVar(...
        shapeByPeriodsToVarsLinIdx ) = 1;

end

% make shape functions at cheb points with shapeFnCheb
TnShape = evalChebClenshaw( chebPts2( numPPShape ), eye( numPPShape  ) );
shapeFnChebCoeffs = TnShape \ eye( numPPShape );

shapeByPeriodsToChebByPeriods = zeros( numPPShape * numPPPeriods, numPPShape * numPPPeriods );
% chebByPeriodsToShapeByPeriods = zeros( numPPShape * numPPPeriods, numPPShape * numPPPeriods );
for nn = 1:numPPPeriods

    numChebNumPIdx = numPPShape * ( nn-1 ) + ( 1 : numPPShape );
    numShapeNumPIdx = numPPShape * ( nn-1 ) + ( 1 : numPPShape );

    shapeByPeriodsToChebByPeriods( numChebNumPIdx, numShapeNumPIdx ) = ...
        shapeFnChebCoeffs;

    % chebByPeriodsToShapeByPeriods( numShapeNumPIdx, numChebNumPIdx ) = TnShape;
end

% evaluate cheb polynomials at points of the time intergration
chebByPeriodsToTimepoints = zeros( numTimePoints, numPPShape * numPPPeriods );
for nn = 1:numPPPeriods
    timePointsPIdx = PPIdxs( nn, 1 ) : PPIdxs( nn, 2 );
    numChebNumPIdx = numPPShape * ( nn-1 ) + ( 1 : numPPShape );

    timePointsChebDomain = chebFwdMap( tvec( timePointsPIdx ) , PPStartStop( nn, : ) );
    Tnp = evalChebClenshaw( timePointsChebDomain, eye( numPPShape  ) );
    chebByPeriodsToTimepoints( timePointsPIdx, numChebNumPIdx ) = ...
        Tnp;
end

varsToChebByPeriods = shapeByPeriodsToChebByPeriods * varsToShapeByPeriods;
varsToTimepoints = chebByPeriodsToTimepoints * varsToChebByPeriods;

% chebByPeriodsToVars = shapeByPeriodsToVar * chebByPeriodsToShapeByPeriods  ;

end