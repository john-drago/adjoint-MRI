function [ c, ceq, gradc, gradceq ] = nonlconSpokesCombine( pSc, opt, spokes, nonlconIneqFuncsSpokes, nonlconEqFuncsSpokes )
% This function will take all of the individually defined nonlcon functions
% and combine them into one function that can be called.

% Assume pSc is scaled to be between -1 and 1

c = [];
ceq = [];
gradc = [];
gradceq = [];

% iterate across nonlinear inequality constraints
nlconIneqFuncs = nonlconIneqFuncsSpokes;
if ~isempty( nlconIneqFuncs )
    for cc = 1:length( nlconIneqFuncs )
        nlconIneqFunc = nlconIneqFuncs{ cc };
        if nargout > 2
            [ ci, gradci ] = nlconIneqFunc( pSc, opt, spokes );

            gradc = [ gradc, gradci ]; %#ok
        else
            ci = nlconIneqFunc( pSc, opt, spokes );
        end
        c = [ c; ci ]; %#ok
        
    end
else
    c = [];
    gradc = [];
end

% iterate across nonlinear equality constraints
nlconEqFuncs = nonlconEqFuncsSpokes;
if ~isempty( nlconEqFuncs )
    for cc = 1:length( nlconEqFuncs )
        nlconEqFunc = nlconEqFuncs{ cc };

        if nargout > 3
            [ ceqi, gradceqi ] = nlconEqFunc( pSc, opt, spokes );
            gradceq = [ gradceq, gradceqi ]; %#ok
        else
            [ ceqi ] = nlconEqFunc( pSc, opt, spokes );
        end

        ceq = [ ceq; ceqi ]; %#ok

    end
else
    ceq = [];
    gradceq = [];
end

end