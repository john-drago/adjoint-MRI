function [ nlconSt ] = nlconAddList( nlconSt, nlconFunction, nlconFunctionOutNum )

if ~isempty( nlconSt.nlconFuncNames )
    isprev = contains( string( nlconSt.nlconFuncNames ), func2str( nlconFunction ) );
else
    isprev = false;
end

if ~isprev
    nlconSt.nlconFuncidx = nlconSt.nlconFuncidx + 1;

    nlconSt.nlconFuncs{ nlconSt.nlconFuncidx, 1 } = nlconFunction;
    nlconSt.nlconFuncNames{ nlconSt.nlconFuncidx, 1 } =...
        func2str( nlconSt.nlconFuncs{ nlconSt.nlconFuncidx } );
    nlconSt.nlconAmts( nlconSt.nlconFuncidx, 1 ) = nlconFunctionOutNum;
end

end