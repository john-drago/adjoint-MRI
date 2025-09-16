function dkydt = dkydt_fn( t, s )
t = t( : ).';

kr = kr_fn( t, s );
ktheta = ktheta_fn( t, s );
kphi = kphi_fn( t, s );
sinktheta = sin( ktheta );
% cosktheta = cos( ktheta );
sinkphi = sin( kphi );
coskphi = cos( kphi );
dkrdt = dkrdt_fn( t, s );
dkthetadt = dkthetadt_fn( t, s );
dkphidt = dkphidt_fn( t, s );
dsinkthetadt = dkthetadt .* cos( ktheta );
% dcoskthetadt = -dkthetadt .* sin( ktheta );
dsinkphidt = dkphidt .* coskphi;
% dcoskphidt = -dkphidt .* sinkphi;

dkydt = ...
    dkrdt .* sinktheta .* sinkphi + ...
    kr .* dsinkthetadt .* sinkphi + ...
    kr .* sinktheta .* dsinkphidt;

end