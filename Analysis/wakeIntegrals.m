function tpintegrals = wakeIntegrals(u_V,runtable,wing_geom)
% performs integration of wake values for WT tests in 3x2 tunnel
% integration is performed down each slice then across the span

sz = size(u_V); % sz(1) = normal, sz(2) = span

% pull out wing geometry
S = wing_geom.b_wing * wing_geom.c_wing;

% integration geom ds,dn
ds = convlength(runtable.dspitot,'in','m');
dn = convlength(runtable.dzpitot,'in','m');
s  = convlength(24,'in','m');
% extract table values

rho = runtable.rho;
Vinf = runtable.Vinf;
qinf = .5*rho*Vinf^2;

dmdot = zeros(1,sz(2)-1);
dJ = zeros(1,sz(2)-1);
dE = zeros(1,sz(2)-1);
Fx = zeros(1,sz(2)-1);
KE = zeros(1,sz(2)-1);
% integration in two steps

% vertical integration
for slice = 1:(sz(2)-1)
    u_Vs = mean([u_V(:,slice) u_V(:,slice+1)],2);
    
    dmdot(slice) = rho*Vinf*trapz(u_Vs-1)*dn;
    dJ(slice) = rho*Vinf^2*trapz(u_Vs.^2-1)*dn;
    dE(slice) = 0.5 * rho*Vinf^3*trapz(u_Vs.^3-1)*dn;
    Fx(slice) =  Vinf*dmdot(slice)- dJ(slice);
    KE(slice) = dE(slice)-.5*Vinf^2*dmdot(slice);
end

% spanwise integration
dmdot_int = sum(dmdot)*ds * s/(sz(2)*ds);
dJ_int = sum(dJ)*ds * s/(sz(2)*ds);
dE_int = sum(dE)*ds * s/(sz(2)*ds);
Fx_int = sum(Fx)*ds * s/(sz(2)*ds);
KE_int = sum(KE)*ds * s/(sz(2)*ds);


% average, integrate in one step
u_Va = mean(u_V,2);

dmdot2D = rho*Vinf*sum(u_Va-1)*dn;
dJ2D    = rho*Vinf^2*sum(u_Va.^2-1)*dn;
dE2D    = 0.5 * rho*Vinf^3*sum(u_Va.^3-1)*dn;
Fx2D    = Vinf*dmdot2D - dJ2D;
KE2D    = -(dE2D-.5*Vinf^2*dmdot2D);

dmdota   = 24 * dmdot2D * ds;
dJa      = 24 * dJ2D * ds ;
dEa      = 24 * dE2D * ds ;
Fxa      = 24 * Fx2D * ds ;
KEa      = 24 * KE2D * ds ;


tpintegrals = table();

% slices
tpintegrals.dmdot = dmdot;
tpintegrals.dJ = dJ;
tpintegrals.dE = dE;
tpintegrals.Fx = Fx;
tpintegrals.KE = KE;
% integrated values
tpintegrals.dmdot_int = dmdot;
tpintegrals.dJ_int = dJ_int;
tpintegrals.dE_int = dE_int;
tpintegrals.Fx_int = Fx_int;
tpintegrals.KE_int = KE_int;

tpintegrals.dCQ_int = dmdot_int/(rho*Vinf*S);
tpintegrals.dCJ_int = dJ_int/(qinf*S);
tpintegrals.dCE_int = dE_int/(qinf*Vinf*S);
tpintegrals.Cx_int  = Fx_int/(qinf*S);

% quasi-2D average
tpintegrals.dmdota = dmdota;
tpintegrals.dJa = dJa;
tpintegrals.dEa = dEa;
tpintegrals.Fxa = Fxa;
tpintegrals.KEa = KEa;

tpintegrals.dCQa = dmdota/(rho*Vinf*S);
tpintegrals.dCJa = dJa/(qinf*S);
tpintegrals.dCEa = dEa/(qinf*Vinf*S);
tpintegrals.Cxa   = Fxa/(qinf*S);


end