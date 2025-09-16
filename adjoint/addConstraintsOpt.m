function opt = addConstraintsOpt( opt, AbSt, AbeqSt, nlconIneqSt, nlconEqSt )

% Deal with possibility that there are no Ab constraints
if AbSt.Abidx == 0
    AbSt.AbConstraintNames = [];
    AbSt.AbConstraintAmts = [];
end
AbSt.AbConstraintCumSum = cumsum( AbSt.AbConstraintAmts );

if AbeqSt.Abidx == 0
    AbeqSt.AbConstraintNames = [];
    AbeqSt.AbConstraintAmts = [];
end
AbeqSt.AbConstraintCumSum = cumsum( AbeqSt.AbConstraintAmts );

% Deal with possibility that there are no nonlinear constraints
if nlconEqSt.nlconFuncidx == 0
    nlconEqSt.nlconFuncs = [];
    nlconEqSt.nlconFuncNames = [];
    nlconEqSt.nlconAmts = [];
end
nlconEqSt.nlconCumSum = cumsum( nlconEqSt.nlconAmts );

if nlconIneqSt.nlconFuncidx == 0
    nlconIneqSt.nlconFuncs = [];
    nlconIneqSt.nlconFuncNames = [];
    nlconIneqSt.nlconAmts = [];
end
nlconIneqSt.nlconCumSum = cumsum( nlconIneqSt.nlconAmts );

% opt.A = full( AbSt.A ); % SQP and GA don't want sparse matrices
opt.A = AbSt.A;
opt.AbConstraintNames = AbSt.AbConstraintNames;
opt.AbConstraintAmts = AbSt.AbConstraintAmts;
opt.AbConstraintCumSum = AbSt.AbConstraintCumSum;
opt.b = AbSt.b;

% opt.Aeq = full( AbeqSt.A ); % SQP and GA don't want sparse matrices
opt.Aeq = AbeqSt.A;
opt.AbeqConstraintNames = AbeqSt.AbConstraintNames;
opt.AbeqConstraintAmts = AbeqSt.AbConstraintAmts;
opt.AbeqConstraintCumSum = AbeqSt.AbConstraintCumSum;
opt.beq = AbeqSt.b;

opt.nlconIneqFuncs = nlconIneqSt.nlconFuncs;
opt.nlconIneqFuncNames = nlconIneqSt.nlconFuncNames;
opt.nlconIneqAmts = nlconIneqSt.nlconAmts;
opt.nlconIneqCumSum = nlconIneqSt.nlconCumSum;

opt.nlconEqFuncs = nlconEqSt.nlconFuncs;
opt.nlconEqFuncNames = nlconEqSt.nlconFuncNames;
opt.nlconEqAmts = nlconEqSt.nlconAmts;
opt.nlconEqCumSum = nlconEqSt.nlconCumSum;

end