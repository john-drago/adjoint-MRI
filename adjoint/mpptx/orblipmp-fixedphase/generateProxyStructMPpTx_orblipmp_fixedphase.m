function popt = generateProxyStructMPpTx_orblipmp_fixedphase( opt, oc )
% This function will generate a proxy opt struct, popt, that will allow a
% proxy optimization problem to be solved.

popt = struct;

%% Define constants
popt.convertMBackToLarmor = opt.convertMBackToLarmor;
popt.convertMtargAwayLarmor = opt.convertMtargAwayLarmor;
popt.constantRotatingFrame = opt.constantRotatingFrame;

popt.fz_mpptx = opt.fz_mpptx;
popt.wz_mpptx = opt.wz_mpptx;
popt.dwxy_mpptx = opt.dwxy_mpptx;
popt.dfxy_mpptx = opt.dfxy_mpptx;
popt.dwxyvec = opt.dwxyvec;

popt.Z0 = opt.Z0;

popt.gyro = opt.gyro;
popt.numZCoils = opt.numZCoils;
popt.numXYCoils = opt.numXYCoils;
popt.numSubj = opt.numSubj;

%% Define control
popt.useGPU = false;

%% Define field
popt.M0 = opt.M0;
popt.Mtarg = opt.Mtarg;
popt.bzsens = opt.bzsens;
popt.b1preal = opt.b1preal;
popt.b1pimag = opt.b1pimag;
popt.db0 = opt.db0;
popt.pos = opt.pos;
popt.zi = opt.zi;
popt.idx = opt.idx;
popt.optroi = opt.optroi;

popt.numPos = opt.numPos;
popt.numPosSubj = opt.numPosSubj;
popt.cum_numPosSubj = opt.cum_numPosSubj;

%% Define waveform
popt.ORSP_i = opt.ORSP_i;
popt.ORSP_f = opt.ORSP_f;
popt.num_ORSP = opt.num_ORSP;

popt.Blip_i = opt.Blip_i;
popt.Blip_f = opt.Blip_f;
if opt.numZCoils > 0
    popt.Blip_Shim_Slew_i = opt.Blip_Shim_Slew_i;
    popt.Blip_Shim_Slew_f = opt.Blip_Shim_Slew_f;
end
popt.Blip_Grad_Slew_i = opt.Blip_Grad_Slew_i;
popt.Blip_Grad_Slew_f = opt.Blip_Grad_Slew_f;

popt.MPSP_i = opt.MPSP_i;
popt.MPSP_f = opt.MPSP_f;
popt.num_MPSP = opt.num_MPSP;
if opt.numZCoils > 0
    popt.MPSP_Shim_Slew_i = opt.MPSP_Shim_Slew_i;
    popt.MPSP_Shim_Slew_f = opt.MPSP_Shim_Slew_f;
end
popt.MPSP_Grad_Slew_i = opt.MPSP_Grad_Slew_i;
popt.MPSP_Grad_Slew_f = opt.MPSP_Grad_Slew_f;

popt.numTimePoints = opt.numTimePoints;
popt.pulseLength = opt.pulseLength;
popt.tvec = opt.tvec;
popt.dt = opt.dt;
popt.dtvec = opt.dtvec;

popt.tORSP = opt.tORSP;
popt.tBlip = opt.tBlip;
popt.tMPSP = opt.tMPSP;

popt.tStORSP = opt.tStORSP;
popt.tEndORSP = opt.tEndORSP;
popt.tStBlip = opt.tStBlip;
popt.tEndBlip = opt.tEndBlip;
popt.tStMPSP = opt.tStMPSP;
popt.tEndMPSP = opt.tEndMPSP;

if opt.numZCoils > 0
    popt.BlipShimSlewTime = opt.BlipShimSlewTime;
    popt.MPSPShimSlewTime = opt.MPSPShimSlewTime;
end
popt.BlipGradSlewTime = opt.BlipGradSlewTime;
popt.MPSPGradSlewTime = opt.MPSPGradSlewTime;

%% Define variables
varArray = {...
    "breal-ORSP", popt.numXYCoils;...
    "bimag-ORSP", popt.numXYCoils;...
    "breal-MPSP", popt.numXYCoils;...
    "bimag-MPSP", popt.numXYCoils;...
    "Gx-Blip", 1;...
    "Gy-Blip", 1;...
    "Gz-Blip", 1;...
    "Gxreal-MPSP", ( 1 );...
    "Gyreal-MPSP", ( 1 );...
    "Gzreal-MPSP", ( 1 );...
    "Gximag-MPSP", ( 1 );...
    "Gyimag-MPSP", ( 1 );...
    "Gzimag-MPSP", ( 1 );...
    "shimmag-MPSP", (popt.numZCoils);...
    "shimph-MPSP", 1;...
    "shimblipscale", 1;...
    };

varNames = string( varArray( :, 1 ) );
if opt.numZCoils == 0
    shimVars = contains(varNames, "shim", 'ignorecase', true );
    varArray = varArray( ~shimVars, : );
end

varNumDictionary = dictionary();

varCumSum = cumsum( cell2mat( [ 1; varArray( :, 2 ) ] ) );
varIdxs = cell( (size( varCumSum, 1 ) - 1), 1 );
for ii = 1 : size( varIdxs, 1 )
    varIdxs{ ii } = uint32( varCumSum( ii ) : ( varCumSum( ii+1 ) - 1 ) ).';
    varNumDictionary( varArray{ ii, 1 } ) = ii;
end
varCumSum = uint32( varCumSum( 1 : (end-1) ) );

varArray = [ varArray, varIdxs ];
varNames = string( varArray( :, 1 ) );
varAmts = uint32( cell2mat( varArray( :, 2 ) ) );

numVars = varCumSum( end ) + varAmts( end ) - 1;

varCumSum = cumsum( varAmts );

popt.varNames = varNames;
popt.varAmts = varAmts;
popt.varNumDictionary = varNumDictionary;
popt.varCumSum = varCumSum;
popt.varIdxs = varIdxs;
popt.numVars = numVars;

%% Assign indices to opt struct for faster indexing in functions
for vv = 1:length( varNames )
    popt.( sprintf( "%s_idx", replace(varNames( vv ), "-", "_" ) ) ) = varIdxs{ vv };
end

%% Define scaling, bounds, and constraints
popt.A = [];
popt.b = [];
popt.Aeq = [];
popt.beq = [];
popt.nonlcon = [];
popt.nlconIneqFuncs = [];
popt.nlconEqFuncs = [];

scVec = getVecFromOpt( 'scVec', popt, opt );
lb = getVecFromOpt( 'lb', popt, opt  );
ub = getVecFromOpt( 'ub', popt, opt  );

popt.scVec = scVec;
popt.lb = lb;
popt.ub = ub;

if isfield( oc, 'popt' )
    if isfield( oc.popt, 'zCoilScale' )
        if opt.numZCoils > 0
            zCoilIdxs = [ popt.shimmag_MPSP_idx ];
            popt.lb( zCoilIdxs ) = oc.popt.zCoilScale * popt.lb( zCoilIdxs );
            popt.ub( zCoilIdxs ) = oc.popt.zCoilScale * popt.ub( zCoilIdxs );
        end
    end
    if isfield( oc.popt, 'rfScale' )
        rfIdxs = [ popt.breal_ORSP_idx; popt.bimag_ORSP_idx; popt.breal_MPSP_idx; popt.bimag_MPSP_idx ];
        popt.lb( rfIdxs ) = oc.popt.rfScale * popt.lb( rfIdxs );
        popt.ub( rfIdxs ) = oc.popt.rfScale * popt.ub( rfIdxs );
    end
end

end

%% Helper functions
% ----------------------------------------------------------------------- %
function vec = getVecFromOpt( vecName, popt, opt )

vec = zeros( popt.numVars, 1 );
vec( popt.breal_ORSP_idx ) = opt.(vecName)( opt.breal_ORSP_idx( 1:opt.num_ORSP:(opt.numXYCoils*opt.num_ORSP) ) );
vec( popt.bimag_ORSP_idx ) = opt.(vecName)( opt.bimag_ORSP_idx( 1:opt.num_ORSP:(opt.numXYCoils*opt.num_ORSP) ) );
vec( popt.breal_MPSP_idx ) = opt.(vecName)( opt.breal_MPSP_idx( 1:opt.num_MPSP:(opt.numXYCoils*opt.num_MPSP) ) );
vec( popt.bimag_MPSP_idx ) = opt.(vecName)( opt.bimag_MPSP_idx( 1:opt.num_MPSP:(opt.numXYCoils*opt.num_MPSP) ) );
vec( popt.Gx_Blip_idx ) = opt.(vecName)( opt.Gx_Blip_idx );
vec( popt.Gy_Blip_idx ) = opt.(vecName)( opt.Gy_Blip_idx );
vec( popt.Gz_Blip_idx ) = opt.(vecName)( opt.Gz_Blip_idx );
vec( popt.Gxreal_MPSP_idx ) = opt.(vecName)( opt.Gxreal_MPSP_idx );
vec( popt.Gyreal_MPSP_idx ) = opt.(vecName)( opt.Gyreal_MPSP_idx );
vec( popt.Gzreal_MPSP_idx ) = opt.(vecName)( opt.Gzreal_MPSP_idx );
vec( popt.Gximag_MPSP_idx ) = opt.(vecName)( opt.Gximag_MPSP_idx );
vec( popt.Gyimag_MPSP_idx ) = opt.(vecName)( opt.Gyimag_MPSP_idx );
vec( popt.Gzimag_MPSP_idx ) = opt.(vecName)( opt.Gzimag_MPSP_idx );
if opt.numZCoils > 0
    vec( popt.shimmag_MPSP_idx ) = opt.(vecName)( opt.shimmag_MPSP_idx );
    vec( popt.shimph_MPSP_idx ) = opt.(vecName)( opt.shimph_MPSP_idx );
    vec( popt.shimblipscale_idx ) = opt.(vecName)( opt.shimblipscale_idx );
end

end
% ----------------------------------------------------------------------- %