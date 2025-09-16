function [ c, ceq, gradc, gradceq ] = nonlconCombine( pSc, opt )
% This function will take all of the individually defined nonlcon functions
% and combine them into one function that can be called.

% Assume pSc is scaled to be between -1 and 1

c = [];
ceq = [];
gradc = [];
gradceq = [];

if nargout > 2
    gradTol = 1e-10;
end

% iterate across nonlinear inequality constraints
nlconIneqFuncs = opt.nlconIneqFuncs;
if ~isempty( nlconIneqFuncs )
    for cc = 1:length( nlconIneqFuncs )
        nlconIneqFunc = nlconIneqFuncs{ cc };
        if nargout > 2
            [ ci, gradci ] = nlconIneqFunc( pSc, opt );

            gradci( abs( gradci ) < gradTol ) = 0;

            gradc = [ gradc, gradci ]; %#ok
        else
            ci = nlconIneqFunc( pSc, opt );
        end
        c = [ c; ci ]; %#ok
        
    end
else
    c = [];
    gradc = [];
end

% iterate across nonlinear equality constraints
nlconEqFuncs = opt.nlconEqFuncs;
if ~isempty( nlconEqFuncs )
    for cc = 1:length( nlconEqFuncs )
        nlconEqFunc = nlconEqFuncs{ cc };

        if nargout > 3
            [ ceqi, gradceqi ] = nlconEqFunc( pSc, opt );
            
            gradceqi( abs( gradceqi ) < gradTol ) = 0;

            gradceq = [ gradceq, gradceqi ]; %#ok
        else
            [ ceqi ] = nlconEqFunc( pSc, opt );
        end

        ceq = [ ceq; ceqi ]; %#ok

    end
else
    ceq = [];
    gradceq = [];
end

end