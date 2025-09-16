function [ AbSt, AbeqSt, nlconIneqSt, nlconEqSt ] = initializeConstraintTrackers()

A = [];
b = [];
Aeq = [];
beq = [];

AbSt = struct;
AbSt.A = A;
AbSt.b = b;
AbSt.AbConstraintNames = strings( 0, 0 );
AbSt.Abidx = 0;

AbeqSt = struct;
AbeqSt.A = Aeq;
AbeqSt.b = beq;
AbeqSt.AbConstraintNames = strings( 0, 0 );
AbeqSt.Abidx = 0;

nlconIneqSt = struct;
nlconIneqSt.nlconAmts = [];
nlconIneqSt.nlconFuncs = {};
nlconIneqSt.nlconFuncNames = strings( 0, 0 );
nlconIneqSt.nlconFuncidx = 0;

nlconEqSt = struct;
nlconEqSt.nlconAmts = [];
nlconEqSt.nlconFuncs = {};
nlconEqSt.nlconFuncNames = strings( 0, 0 );
nlconEqSt.nlconFuncidx = 0;

end