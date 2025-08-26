function [ constBin, maxConstrViol ] = runConstraintsConvergence( pSc, const, constrTol )

if nargin < 3
    constrTol = 5.0e-2;
end

maxConstrViol = -inf;

if ~isempty( const.A )
    cA = const.A * pSc - const.b;
    cA_bin = all( cA < constrTol );
    boundsBin = cA_bin;
    maxConstrViol = max( [ cA; maxConstrViol ] );
else
    boundsBin = true;
end

if boundsBin
    if ~isempty( const.Aeq )
        cAeq = const.Aeq * pSc - const.beq;
        cAeq_bin = all( ( cAeq < constrTol ) & ( cAeq > ( -constrTol ) ) );
        boundsBin = cAeq_bin;
        maxConstrViol = max( [ abs( cAeq ); maxConstrViol ] );
    else
        boundsBin = true;
    end
else
    constBin = false;
    return;
end

if boundsBin
    [ cnlconIneq, cnlconEq ] = const.nonlcon( pSc );

    if ~isempty( cnlconIneq )
        cnlconIneq_bin = all( cnlconIneq < constrTol );
        boundsBin = cnlconIneq_bin;
        maxConstrViol = max( [ cnlconIneq; maxConstrViol ] );
    else
        boundsBin = true;
    end

    if boundsBin
        if ~isempty( cnlconEq )
            cnlconEq_bin = all( ( cnlconEq < constrTol ) & ( cnlconEq > ( -constrTol ) ) );
            boundsBin = cnlconEq_bin;
            maxConstrViol = max( [ abs( cnlconEq ); maxConstrViol ] );
        else
            boundsBin = true;
        end
    else
        constBin = false;
        return;
    end
else
    constBin = false;
    return;
end

if boundsBin
    constBin = true;
    return;
else
    constBin = false;
    return;
end

end