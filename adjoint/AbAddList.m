function AbSt = AbAddList( AbSt, constraintName, Anew, bnew )

if ~isempty( AbSt.AbConstraintNames )
    isprev = contains( string( AbSt.AbConstraintNames ), constraintName );
else
    isprev = false;
end

if ~isprev
    AbSt.Abidx = AbSt.Abidx + 1;

    AbSt.A = [ AbSt.A; Anew ];
    AbSt.b = [ AbSt.b; bnew ];

    AbSt.AbConstraintNames( AbSt.Abidx, 1 ) = constraintName;
    AbSt.AbConstraintAmts( AbSt.Abidx, 1 ) = uint32( size( Anew, 1 ) );
end

end