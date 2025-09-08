function Mvec = blochSimVec( simSt, tol )
% blochSimVec
% Vectorized Bloch-equation simulation using MATLAB’s ODE solver in the
% rotating frame. Computes magnetization dynamics for a set of spatial
% positions given gradient, RF, shim, and sensitivity inputs.
%
% INPUTS:
%   simSt : struct with fields RF(Nb1×T, complex), Grad(3×T), Shim(Nshim×T),
%           gyro (scalar), tvec(1×T), pos(N×3), M0(N×3), 
%           b1preal(N×Nb1), b1pimag(N×Nb1), bzsens(N×Nz), db0(N×1)
%   tol   : ODE rel/abs tolerance (default 1e-6)
%
% OUTPUT:
%   Mvec  : final magnetization state matrix (N×3)

arguments
    simSt
    tol = 1e-6;
end

dMdtSt = struct;
dMdtSt.RF = simSt.RF;
dMdtSt.Grad = simSt.Grad;
dMdtSt.gyro = simSt.gyro;
dMdtSt.tvec = simSt.tvec;
dMdtSt.xyz = simSt.pos;
dMdtSt.M0 = simSt.M0;
dMdtSt.Shim = simSt.Shim;
dMdtSt.BxySensCoil = simSt.b1preal + 1j * simSt.b1pimag;
dMdtSt.BzSensCoilz = simSt.bzsens;
dMdtSt.DB0Vec = simSt.db0;

tend = simSt.tvec(end);

dMdtSt.M0xyzstack = reshape( dMdtSt.M0.', [ numel( dMdtSt.M0 ), 1 ] );

%% Sim Loop
options = odeset('RelTol', tol, 'AbsTol', tol);

%% Sim Loop
if tend == 0
    Msim = M0xyzstack;
else
    [~, Msim] = ode45(...
        @(t,M) dMdtMP(t,M,dMdtSt),...
        [0, tend/2, tend], dMdtSt.M0xyzstack, options);

    Msim = Msim(end,:).';
end

Mvec = reshape( Msim, [ size( dMdtSt.M0, 2 ), size( dMdtSt.M0, 1 ) ] ).';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dMdt = dMdtMP(t, Mt, dMdtSt)
% in the rotating frame at Larmor frequency

% This is how you would calculate the fields from the B1z coils in the
% rotating frame. As there is a massive difference between the Larmor
% Frequency at 7 T and the frequency of our z-directed field, this will
% take forever to compute.

% Z Fields
Shim = interp1( dMdtSt.tvec, dMdtSt.Shim.', t ).';
BzShim = dMdtSt.BzSensCoilz * ( Shim );
Grad = interp1( dMdtSt.tvec, dMdtSt.Grad.', t ).';
BzGrad = dMdtSt.xyz * Grad;
DB0 = dMdtSt.DB0Vec;

% XY Fields
RF = interp1( dMdtSt.tvec, dMdtSt.RF.', t).';
BxyRF = dMdtSt.BxySensCoil * RF;

Bt = [...
    real(BxyRF),...
    imag(BxyRF),...
    (BzGrad + BzShim + DB0)].';

Mt = reshape(Mt, [3, size(Bt,2)]);

dMdt = dMdtSt.gyro * (cross(Mt, Bt));

dMdt = reshape(dMdt, [3*size(dMdt,2),1]);

if any(isnan(dMdt))
    error('Returned NaN value in Bloch Sim');
end
end