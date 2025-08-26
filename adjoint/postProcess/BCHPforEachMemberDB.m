function [ BCmagVec, FABCHP, MBCHP, BCHP_pulsePower ] = BCHPforEachMemberDB( val )
% This function will create birdcage hardpulse for each member in the
% database

% Initialize return structs
FABCHP = zeros( val.cum_numPosSubj( end ), 1  );
MBCHP = zeros( val.cum_numPosSubj( end ), 3 );
BCmagVec = zeros( val.numSubj, 1 );
BCHP_pulsePower = zeros( val.numSubj, 1 );

for nn = 1:val.numSubj
    if nn == 1
        NonOptIdxs = 1:( val.cum_numPosSubj(nn));
    else
        NonOptIdxs = ( val.cum_numPosSubj(nn-1)+1 ):( val.cum_numPosSubj(nn) );
    end
    
    % Calculate Mxy in brain that arises from one volt
    optIdxs = NonOptIdxs( val.optroi(NonOptIdxs) );
    bSensMagBrain = sqrt( abs( val.b1p( optIdxs, : ) ).^2 + val.db0(optIdxs).^2 );
    bSensVecBrain = [ real( val.b1p( optIdxs, : ) ), imag( val.b1p( optIdxs, : ) ), val.db0(optIdxs) ];
    bSensUnitBrain = bSensVecBrain ./ bSensMagBrain;
    bSensUnitBrain( isnan(bSensUnitBrain) ) = 0;
    thetaUnitBrain = -val.gyro * val.tBCHP * bSensMagBrain;
    M0VecBrain = zeros( size( bSensVecBrain ) );
    M0VecBrain( :, 3 ) = 1;
    MsensBrain = calcRRF_vec( M0VecBrain, bSensUnitBrain,  thetaUnitBrain );
    MxysensBrain = MsensBrain( :, 1 ) + 1j * MsensBrain( :, 2 );
    Mtarg = val.Mtarg( optIdxs, : );
    targFA = real( acosd( Mtarg( :, 3 ) ) );
    
    % scale birdcage voltage to get target flip angle
    initBCmag = mean( targFA ./ real( asind( abs( MxysensBrain ) ) ) );

    [BCmag, ~] = optBCmagBCHP( initBCmag, targFA, bSensVecBrain, M0VecBrain, val.gyro, val.tBCHP );

    % Calculate Mxy in body that arises from one volt
    bSensVecBody = [ real( val.b1p( NonOptIdxs, : ) ), imag( val.b1p( NonOptIdxs, : ) ), val.db0(NonOptIdxs) ];
    M0VecBody = zeros( size( bSensVecBody ) );
    M0VecBody( :, 3 ) = 1;

    [ ~, FABCHP_subj, MBCHP_subj ] = getMeanFlipAngle( BCmag, bSensVecBody, M0VecBody, val.gyro, val.tBCHP );

    MBCHP( NonOptIdxs, : ) = MBCHP_subj;
    FABCHP( NonOptIdxs ) = FABCHP_subj;

    BCmagVec( nn ) = BCmag;
    
    BCHPwaveform = BCmag * ones( 1,length(val.tvec) );
    tend = val.tvec( end );
    dt = val.dt;
    Z0 = val.Z0;
    BCHP_pulsePower( nn ) = ...
        1/( tend * Z0 ) * sum( dt * BCHPwaveform.^2 );

end

end