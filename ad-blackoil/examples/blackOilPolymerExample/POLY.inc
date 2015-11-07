NOECHO
-- Polymer properties (keywords starting with PLY)
-- PLYVISC: Viscosity multiplier. Defined as PLYVISC=u/uw, where uw is water viscosity. uw@50 C is 0.5 mPa*s (cP)
--          Tables defined twice as NTSFUN=2
PLYVISC
0    1.0
0.5  3.5
1.0  6.4
1.5  12.4
2.0  24.7
3.0  50.0 /


-- PLYROCK gives fraction of inaccessible pore volume (IPV, must be less than Swi)
-- residual resistivity factor (RRF), density of rock incl. PV (P-sandgrain=2.7 g/cm3 + poro=0.185)
-- desorpsjon av polymer (paa=1,av=2), og RRF. ved maksimalt mengde adsorbert polymer
-- IPV: Eclipse krever at IPV<Swcr. Vi har derfor satt denne til 0.199. Our recommended base-case value is 0.25
-- RRF: 33% reduction in water. perm (33% was observed for Oerrell)
-- AI: adsorbtion index: 1=desorption may occur, 2= NO desorption
-- Max. adsorbert polymer (kg/kg),

PLYROCK
--IPV  RRF  dens  AI  max ads
 0.05 1.3  2650   2 0.000020 /


--PLYADS: Adsorption isotherm
PLYADS
-- conc  adsorb-conc
 0.000   0.000000
 0.250   0.000011
 0.500   0.000015
 0.750   0.000017
 1.000   0.000018
 1.250   0.000019
 1.500   0.000019
 1.750   0.000020
 2.000   0.000020
 3.000   0.000020 /

--Polymer: Todd-Longstaff mixing parameter. 1 = full mixing
PLMIXPAR
 1.0  /

--Max. concentration to be used in the mixing calculation (set equal to pol. inj. conc.)
PLYMAX
 3.00  0  /

-- PLYSHEAR: defined by a viscosity reduction factor (VRF) as a function of water phase velocity (Vw).
-- Omit shear as CPU-time is significantly reduced and shear thinning is only significant for injectivity which
--  will have to be studied by other tools than Eclipse
PLYSHEAR
-- ---- Vw.            VRF
0           1
0.0000216   0.9
0.0000432   0.85
0.0000864   0.8
0.0001296   0.79
0.0001728   0.78
0.000216    0.77
0.000432    0.74
0.000648    0.71
0.000864    0.69
0.00108     0.68
0.001296    0.67584
0.001512    0.6688
0.001728    0.6688
0.001944    0.66176
0.00216     0.66176
0.00432     0.6336
0.00648     0.61952
0.00864     0.61248
0.0108      0.5984
0.01296     0.59136
0.01512     0.58432
0.01728     0.57728
0.01944     0.57728
0.0216      0.57024
0.216       0.47872
2.16        0.40832   /

ECHO
