function [ opt, snopt ] = prepareSNOPTFunctions( opt, snopt )

%% Define objective function
opt.costFnOpt = @( pSc ) objectiveFunction( pSc, opt );

%% Define nonlcon function
opt.nonlconSNOPT = @( pSc ) nlconFunction( pSc, opt.nonlcon );

end

%% Helper Functions
% ----------------------------------------------------------------------- %
function [ c, ceq, gradc, gradceq ] = nlconFunction( pSc, nlcon )

if nargout > 2
    [ c, ceq, dc, dceq ] = nlcon( pSc );
    gradc = transpose( dc );
    gradceq = transpose( dceq );
    
    gradTol = 1e-10;
    gradc( abs( gradc ) < gradTol ) = 0;
    gradceq( abs( gradceq ) < gradTol ) = 0;

else
    [ c, ceq ] = nlcon( pSc );
end

end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function [ J, dJdp ] = objectiveFunction( pSc, opt )

costFn = opt.runCostFunction;

if nargout > 1
    [ J, dJdp ] = costFn( pSc( : ), opt );
    
    % dJdp = transpose( gradpSc_J );
    % dJdp = gradpSc_J;

    gradTol = 1e-10;
    dJdp( abs( dJdp ) < gradTol ) = 0;

else
    J = costFn( pSc( : ), opt );
end

end
% ----------------------------------------------------------------------- %