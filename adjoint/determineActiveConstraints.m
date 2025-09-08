function [ AbActive, AbeqActive, lbubActive, nonlconIneqActive, nonlconEqActive ] = ...
    determineActiveConstraints( pSc, opt, tolConstraint )
% This function will determine which constraints are active

%% A, b constraints
if ~isempty( opt.A )
    AbVals = ( opt.A * pSc - opt.b );
    activeAbidxs = find( AbVals > ( tolConstraint ) );

    if ~isempty( activeAbidxs )

        AbConstraintStart = [ 1; ( opt.AbConstraintCumSum( 1:( end - 1) ) + 1 ) ];
        AbConstraintFinal = opt.AbConstraintCumSum;

        AbActive = cell( length( activeAbidxs ), 4 );

        for aa = 1:length( activeAbidxs )
            constLogical = ( activeAbidxs( aa ) >= AbConstraintStart ) & ( activeAbidxs( aa ) <= AbConstraintFinal );

            AbActive{ aa, 1 } = opt.AbConstraintNames( constLogical, 1 );
            AbActive{ aa, 2 } = uint32( activeAbidxs( aa ) - AbConstraintStart( constLogical ) + 1 );
            AbActive{ aa, 3 } = uint32( activeAbidxs( aa ) );
            AbActive{ aa, 4 } = AbVals( activeAbidxs( aa ) );
        end

    else
        AbActive = [];
    end
else
    AbActive = [];
end

%% Aeq, beq constraints
if ~isempty( opt.Aeq )
    AbeqVals = abs( opt.Aeq * pSc - opt.beq );
    activeAbeqidxs = find( AbeqVals > ( tolConstraint ) );
    if ~isempty( activeAbeqidxs )

        AbeqConstraintStart = [ 1; ( opt.AbeqConstraintCumSum( 1:( end - 1) ) + 1 ) ];
        AbeqConstraintFinal = opt.AbeqConstraintCumSum;

        AbeqActive = cell( length( activeAbeqidxs ), 3 );

        for aa = 1:length( activeAbeqidxs )
            constLogical = ( activeAbeqidxs( aa ) >= AbeqConstraintStart ) & ( activeAbeqidxs( aa ) <= AbeqConstraintFinal );

            AbeqActive{ aa, 1 } = opt.AbeqConstraintNames( constLogical, 1 );
            AbeqActive{ aa, 2 } = uint32( activeAbeqidxs( aa ) - AbeqConstraintStart( constLogical ) + 1 );
            AbeqActive{ aa, 3 } = uint32( activeAbeqidxs( aa ) );
            AbeqActive{ aa, 4 } = AbeqVals( activeAbeqidxs( aa ) );
        end
    else
        AbeqActive = [];
    end
else
    AbeqActive = [];
end

%% lb, ub constraints
if (~isempty( opt.ub )) || (~isempty( opt.lb ))

    activelbubidxs = find( (abs( pSc - opt.ub ) < ( tolConstraint )) | (abs( pSc - opt.lb ) < ( tolConstraint )) );
    if ~isempty( activelbubidxs )
        lbubConstraintStart = [ 1; ( opt.varCumSum( 1:( end - 1) ) + 1 ) ];
        lbubConstraintFinal = opt.varCumSum;

        lbubActive = cell( length( activelbubidxs ), 3 );

        for aa = 1:length( activelbubidxs )
            constLogical = ( activelbubidxs( aa ) >= lbubConstraintStart ) & ( activelbubidxs( aa ) <= lbubConstraintFinal );

            lbubActive{ aa, 1 } = opt.varNames( constLogical, 1 );
            lbubActive{ aa, 2 } = uint32( activelbubidxs( aa ) - lbubConstraintStart( constLogical ) + 1 );
            lbubActive{ aa, 3 } = uint32( activelbubidxs( aa ) );
            lbubActive{ aa, 4 } = pSc( activelbubidxs( aa ) );
        end
    else
        lbubActive = [];
    end
else
    lbubActive = [];
end

%% nonlcon constraints
if ~isempty( opt.nonlcon )
    [ c, ceq ] = opt.nonlcon( pSc );

    % Deal with inequality constraints
    if ~isempty( c )

        activenonlconIneqidxs = find( c > ( tolConstraint ) );

        if ~isempty( activenonlconIneqidxs )
            nonlconIneqConstraintStart = [ 1; ( opt.nlconIneqCumSum( 1:( end - 1) ) + 1 ) ];
            nonlconIneqConstraintFinal = opt.nlconIneqCumSum;

            nonlconIneqActive = cell( length( activenonlconIneqidxs ), 3 );

            for aa = 1:length( activenonlconIneqidxs )
                constLogical = ( activenonlconIneqidxs( aa ) >= nonlconIneqConstraintStart ) & ( activenonlconIneqidxs( aa ) <= nonlconIneqConstraintFinal );

                nonlconIneqActive{ aa, 1 } = opt.nlconIneqFuncNames{ constLogical, 1 };
                nonlconIneqActive{ aa, 2 } = uint32( activenonlconIneqidxs( aa ) - nonlconIneqConstraintStart( constLogical ) + 1 );
                nonlconIneqActive{ aa, 3 } = uint32( activenonlconIneqidxs( aa ) );
                nonlconIneqActive{ aa, 4 } = c( activenonlconIneqidxs( aa ) );
            end
        else
            nonlconIneqActive = [];
        end
    else
        nonlconIneqActive = [];
    end

    % Deal with equality constraints
    if ~isempty( ceq )
        ceqAbs = abs( ceq );
        activenonlconEqidxs = find( ceqAbs > ( tolConstraint ) );

        if ~isempty( activenonlconEqidxs )
            nonlconEqConstraintStart = [ 1; ( opt.nlconEqCumSum( 1:( end - 1) ) + 1 ) ];
            nonlconEqConstraintFinal = opt.nlconEqCumSum;

            nonlconEqActive = cell( length( activenonlconEqidxs ), 3 );

            for aa = 1:length( activenonlconEqidxs )
                constLogical = ( activenonlconEqidxs( aa ) >= nonlconEqConstraintStart ) & ( activenonlconEqidxs( aa ) <= nonlconEqConstraintFinal );

                nonlconEqActive{ aa, 1 } = opt.nlconFuncNames{ constLogical, 1 };
                nonlconEqActive{ aa, 2 } = uint32( activenonlconEqidxs( aa ) - nonlconEqConstraintStart( constLogical ) + 1 );
                nonlconEqActive{ aa, 3 } = uint32( activenonlconEqidxs( aa ) );
                nonlconEqActive{ aa, 4 } = ceqAbs( activenonlconEqidxs( aa ) );
            end
        else
            nonlconEqActive = [];
        end
    else
        nonlconEqActive = [];
    end

else
    nonlconIneqActive = [];
    nonlconEqActive = [];
end

end