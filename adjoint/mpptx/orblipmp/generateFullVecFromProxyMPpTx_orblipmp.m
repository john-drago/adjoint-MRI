function pSc = generateFullVecFromProxyMPpTx_orblipmp( opt, wvpopt )

% Initialize p struct
p = zeros( opt.numVars, 1 );

% Add ORSP
p( opt.breal_ORSP_idx ) = reshape( wvpopt.brealphasor( :, opt.ORSP_i:opt.ORSP_f ).', [ opt.num_ORSP*opt.numXYCoils, 1 ] );
p( opt.bimag_ORSP_idx ) = reshape( wvpopt.bimagphasor( :, opt.ORSP_i:opt.ORSP_f ).', [ opt.num_ORSP*opt.numXYCoils, 1 ] );
p( opt.breal_MPSP_idx ) = reshape( wvpopt.brealphasor( :, opt.MPSP_i:opt.MPSP_f ).', [ opt.num_MPSP*opt.numXYCoils, 1 ] );
p( opt.bimag_MPSP_idx ) = reshape( wvpopt.bimagphasor( :, opt.MPSP_i:opt.MPSP_f ).', [ opt.num_MPSP*opt.numXYCoils, 1 ] );

% Add Blip
p( opt.Gx_Blip_idx ) = wvpopt.G_Blip( 1 );
p( opt.Gy_Blip_idx ) = wvpopt.G_Blip( 2 );
p( opt.Gz_Blip_idx ) = wvpopt.G_Blip( 3 );

if opt.numZCoils > 0
    p( opt.shim_Blip_idx ) = wvpopt.shim_Blip;
end

% Add MPSP
p( opt.Gxreal_MPSP_idx ) = wvpopt.G_real_MPSP( 1 );
p( opt.Gyreal_MPSP_idx ) = wvpopt.G_real_MPSP( 2 );
p( opt.Gzreal_MPSP_idx ) = wvpopt.G_real_MPSP( 3 );
p( opt.Gximag_MPSP_idx ) = wvpopt.G_imag_MPSP( 1 );
p( opt.Gyimag_MPSP_idx ) = wvpopt.G_imag_MPSP( 2 );
p( opt.Gzimag_MPSP_idx ) = wvpopt.G_imag_MPSP( 3 );

if opt.numZCoils > 0
    p( opt.shimreal_MPSP_idx ) = wvpopt.shimreal_MPSP;
    p( opt.shimimag_MPSP_idx ) = wvpopt.shimimag_MPSP;
end

pSc = p ./ opt.scVec;

end