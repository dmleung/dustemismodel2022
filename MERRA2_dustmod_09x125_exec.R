#############################################################
#############################################################

# dust emission code version 1: 2006 MERRA2 1x1 data for dust emissions
# extract N Africa for calculation only (or specify other spatial domains)
# allow 1-hourly, 2-hourly data or other time resolution simulations (can skip some timestep and multiply dt at the end if one think e.g. 2-hourly sampling is good enough)
# allow showing actual dust emissions and annually scaled emissions (Tg / yr)
# set up spatial and temporal domain

# 28 Oct 2021
# remaining issues:
# 1. Now all MERRA2 variables are 2-hourly and are directly averaged instead of area-weighted averaged. Will fix this next time.
# 2. Check if the boundaries have too much emissions and whether it is needed to multiply by a land fraction to discount those emissions, or else the upscaled variables over the seashores will misbehave (very large or very small).
# 3. Now I have a big time loop for everything; but if later on find that the loop takes too long, one can easily change the script to several smaller time loops and run one by one.
# 4. indices are too bad when extracting MERRA2 variables. Please fix it (8 Oct 2021)

# set working directory
setwd('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR')
# load libraries 
library(ncdf4); library(fields); library(maps); library(Metrics); library(pracma)
# load functions
source("get_geo.R"); source("get_met.R"); source("sptial_plot_fns.R")
source("dust_research_fns.R")
# load color scheme
load("WBGYR_scheme.RData"); load('TEMP_DIFF.RData')

##################
# user input here
# specialize temporal domain as YYYYMMDD
startdate = 20060101   # YYYYMMDD
enddate = 20061231   # YYYYMMDD
# specialize temporal resolution (dt = 1 as hourly, dt = 3 as 3-hourly etc.)
dt = 2   # try 2-hourly
# specialize spatial domain as indices
indi = 1:288; indj = 36:186 
# choose clay dataset: 'FAO' or 'SG' (SoilGrids)
clay_dataset = 'FAO'
# choose the wind scale for dust emission, including default (aerodynamic + nonneutral) ustar (default), aeolian + nonneutral star (aeo+nonneu), and aeolian + neutral ustar (aeo+neu)
wind.scale='default'
# choose land cover averaging method: LC0 (rock and veg) or LC1 (rock, veg and mix)
LC = 'LC0'
# choose the averaging method of the coarse input variables: 'mean' or 'aw_mean' (area-weighted mean)
avg.method = 'aw_mean'
# choose whether default Kok et al. (2014) scheme (K14) or intermittency (K14+Comola et al., 2019) is used. 'FALSE' is K14, 'TRUE' is K14+C19 intermittency scheme.
use.intermittency = FALSE
# fragmentation exponent limit (CESM simulations can blow up when frag_expt is too large when using impact threshold for dust emission equation)
frag_expt_lim = 3   
##################
# calculate model parameters based on user inputs
date_vec = make.date.vec(startdate, enddate)  # make date vec
time = 1:(length(date_vec)*24/dt)  # a vector of timesteps
indt = seq(dt, 24, by=dt)   # hour of day
# load 0.9x1.25 CESM lat lon
nc = nc_open('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR/griddata_0.9x1.25_070212.nc')
#lonc1 = ncvar_get(nc, 'LONW')[,1]
LON.C1 = seq(-179.375,179.375,by=1.25)[indi]; LAT.C1 = ncvar_get(nc, 'LATS')[1,2:192][indj]


###################
# constant parameters
C_tune = 0.05   # global tuning constant
gamma = 5e-4;   # soil aggregation effect on dust emission threshold
u_star_st0 = 0.16;  # m s-1, minimum value of u_star_st
g = 9.81        # m/s^2, gravitational acceleration 
dns_slt = 2650  # kg/m3, soil particle density
D_p = 130e-6    # m, soil median diameter
z0s = D_p / 15    # m, smooth roughness length, White (2006) equation
X = 10            # m, we assume interobstacle distance ~ 10 m
LAI_thr = 1   # Mahowald et al. (2010) set as 0.3 but we changed to 1
U_0 = 0.32  # SSR in the immediate lee of a plant for Okin's scheme, dimensionless
c = 4.8     # e-folding distance SSR recovery for Okin's scheme, dimensionless
k = 0.387     # von Karman constant for high Reynolds number flows
SHR_CONST_RHOFW = 1000   # kg / m3, pure water density
B_it = 0.81 # impact threshold to dry fluid threshold ratio
c_p = 1005   # J/kg/K  specific heat of air at constant pressure
rho_a0 = 1.225; # kg m-3, standard atmospheric density
ustar_st0 = 0.16;  # m s-1, minimum value of standardized fluid threshold
C_d0 = 4.4e-5; C_e = 2.0; C_alpha = 2.7 # dust emission coefficient parameters

###################
# load required time-invariant datasets
.env = new.env()   # set a new environment for storing data temporarily
# load GLCNMO LULC data
load('GLCNMO_frc_area_1x1.RData',envir=.env); ls(.env)
frc.r.C1 = .env$frc.r.C1[indi,indj]; frc.veg.C1 = .env$frc.veg.C1[indi,indj]
frc.mix.C1 = .env$frc.mix.C1[indi,indj]
#A_r = frc.r.C1 / (frc.r.C1+frc.veg.C1+frc.mix.C1)
#A_veg = frc.veg.C1 / (frc.r.C1+frc.veg.C1+frc.mix.C1)
#A_mix = frc.mix.C1 / (frc.r.C1+frc.veg.C1+frc.mix.C1)
#rm('frc.r.C1','frc.veg.C1','frc.mix.C1')
#plot.field(frc.mix.C1, LON.C1, LAT.C1, col=WBGYR)
# load landfilt 
load('landfilt_1x1.RData',envir=.env)
landfilt.C1 = .env$landfilt.C1[indi,indj]
# get Prigent 2005 roughness length and rock drag partition factor
load('Pr05_z0.min_1x1.RData',envir=.env)
z0.r = .env$z0.C1[indi,indj]/100   # in cm, divide by 100 to convert to m
# in this formulation, set globally constant z0s because soil median diameter ~ 130 um
f_eff.r = 1 - ( log(z0.r/z0s) / log(0.7*(X/z0s)^0.8) )  
f_eff.r = f_eff.r*landfilt.C1   # take away values in Caspian sea
#plot.field(f_eff.r, LON.C1, LAT.C1)
# image(f_eff.r)
# get porosity  (issue: MERRA2 porosity may be different from CESM's effective porosity)
load('MERRA2_porosity_1x1.RData', envir=.env)
poros = .env$poros.C1[indi,indj] # porosity (m3 / m3, or volumetric soil water) 
#image(poros)
#plot.field(poros, LON.C1, LAT.C1)
# get bulk density of dry soil material [kg/m^3]
bulk_den = (1 - poros)*dns_slt 
#plot.field(bulk_den, LON.C1, LAT.C1, legend.mar=4)
# use CESM's FAO clay fraction
for (aa in 1) {
  if (clay_dataset == 'FAO') {
    load('f_clay_FAO_1x1.RData',envir=.env)  # CESM's FAO (2012) clay
    f_clay.FAO.1x1 = .env$f_clay.FAO.C1  # 0-1, fraction of #plot.field(f_clay.FAO.05x06, LON.CESM, LAT.CESM, type='def', zlim=c(0,0.65))
    #plot.field(f_clay.FAO.1x1, LON.C1, LAT.C1)
    f_clay.FAO.1x1 = f_clay.FAO.1x1[indi,indj] * landfilt.C1
    #plot.field(f_clay.FAO.1x1, LON.C1, LAT.C1)
    # threshold gravimetric moisture and clay limitation on dust emission
    gwc_thr = 0.01*(0.17*f_clay.FAO.1x1*100 + 0.0014*(f_clay.FAO.1x1*100)^2)*landfilt.C1
    #plot.field(gwc_thr, LON.C1, LAT.C1)
    mss_frc_cly_vld = f_clay.FAO.1x1 * landfilt.C1
    mss_frc_cly_vld[which(f_clay.FAO.1x1>0.2)] = 0.2
    #plot.field(mss_frc_cly_vld, LON.C1, LAT.C1)
    print('used FAO clay data')
  }  else if (clay_dataset == 'SG') {
    load('f_clay_SoilGrids_1x1.RData',envir=.env)
    f_clay.SG.1x1 = .env$f_clay.SG.C1
    f_clay.SG.1x1 = f_clay.SG.1x1[indi,indj]
    #plot.field(f_clay.SG.1x1, LON.C1, LAT.C1, type='def', zlim=c(0,0.65))
    # threshold gravimetric moisture and clay limitation on dust emission
    gwc_thr = 0.01*(0.17*f_clay.SG.1x1*100 + 0.0014*(f_clay.SG.1x1*100)^2)
    #plot.field(gwc_thr, LON.C1, LAT.C1)
    mss_frc_cly_vld = f_clay.SG.1x1
    mss_frc_cly_vld[which(mss_frc_cly_vld>0.2)] = 0.2
    #plot.field(mss_frc_cly_vld, LON.C1, LAT.C1)
    print('used SoilGrids clay data')
  }
}

###################
# scrap files for selected dates
files = NULL   # character strings storing the file names
for (dd in 1:length(date_vec)) {
  date_vec[dd]
  files = c(files, Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/REGRIDDED_MERRA2/',avg.method,'/merra2_regridded_met_1x1_', date_vec[dd], '.RData', sep='')))  
}

###################
# start the main loop
# for days in d
d = 1   # d-th day

for (d in 1:length(date_vec)) {# d-th day
  print(paste('day ', d, ': ', date_vec[d], sep=''))
  #print('extract MERRA2 data', sep='')
  ###################
  # get 3D data in time and space, subsetted as regional data with variable time resolutions
  i1 = indi[1]; j1 = indj[1]; t1 = indt[1]
  il = length(indi); jl = j1+length(indj)-1; tl = t1+length(indt)
  
  #ptm = proc.time()
  #nc = nc_open(files_flx[d])#; names(nc$var)#, start=c(11,26), count=c(5,5)
  load(files[d], envir=.env)#; ls(.env)
  USTAR = .env$USTAR.C1[i1:il,j1:jl,] # m/s
  Z0M   = .env$Z0M.C1[i1:il,j1:jl,]     # m
  RHOA  = .env$RHOA.C1[i1:il,j1:jl,]     # kg/m3
  PBLH  = .env$PBLH.C1[i1:il,j1:jl,]     # kg/m3
  #DISPH = .env$DISPH   # m
  T10M  = .env$T10M.C1[i1:il,j1:jl,]     # K
  #U10   = .env$U10.C1[i1:il,j1:jl,]  # m/s. Define U10 as the total wind speed of U10M and V10M. Don't confuse U10 with U10M. 
  #proc.time() - ptm
  SFMC  = .env$SFMC.C1[i1:il,j1:jl,]
  LAI   = .env$LAI.C1[i1:il,j1:jl,]
  SHLAND = .env$SHLAND.C1[i1:il,j1:jl,]
  #plot.field(DISPH[,,3], LON.C1, LAT.C1)
  #plot.field(RHOA[,,3], LON.C1, LAT.C1)
  #plot.field(Z0M[,,3], LON.C1, LAT.C1)
  #plot.map.log((Z0M[,,3]),LON.C1, LAT.C1,type='def',zlim=c(-4,0), legend.mar=5)
  
  #print('finished loading MERRA2 met fields; now calculate hybrid drag partition and hybrid roughness length')
  
  ###################
  # get drag partition factors and hybrid roughness length (later on, make this as a separate time loop independent of dust emission calculation?)
  
  # get vegetation roughness length following Okin's scheme using LAI
  # first use LAI to get bare land fraction using Natalie Mahowald's equation
  # and vegetation drag partition using Okin-Pierre parameterization
  # Natalie's bare land fraction approximation
  f_bare = 1 - LAI/LAI_thr; f_bare[which(f_bare<0)] = NaN
  f_bare = sweep(f_bare, c(1,2), landfilt.C1, '*')   # optional
  #f_bare_0.3 = 1 - LAI[,,3]/0.3; f_bare_0.3[which(f_bare_0.3<0)] = NaN
  #plot.field(f_bare[,,6], LON.C1, LAT.C1)
  #plot.field(f_bare_0.3, LON.C1, LAT.C1)
  #plot.field.log(LAI[,,6], LON.C1, LAT.C1, zlim=c(-2,1), legend.mar=5)
  
  # make a filter for dust emission sources
  #emisfilt = f_bare/f_bare    # this filter is 3D and changes with time; or one can make a annual, 2D one. 
  #image(emisfilt)
  
  # Okin's vegetation drag partition
  # Using Caroline Pierre (2014)'s formulation of Okin's scheme
  K = pi/2 * (1/LAI - 1); K[which(K<0)] = 0
  SSR = (K + U_0 * c) / (K + c)  # Pierre-Okin shear stress ratio (SSR)
  #plot.field(SSR[,,3], LON.C1, LAT.C1, type='def', zlim=c(0.6,0.9))
  
  # combine rock and vegetation drag partition factors
  # then combine drag partition factors
  #F_eff.LC = array(NaN, dim=c(length(indi), length(indj), length(indt)))
  #F_eff.LC2 = array(NaN, dim=c(length(indi), length(indj), length(indt)))
  #for (h in 1:length(indt)) {
  #F_eff.LC1[,,h] = (frc.r.C1*f_eff.r^3 + frc.veg.C1*SSR[,,h]^3 + frc.mix.C1*(SSR[,,h]*f_eff.r)^3)^(1/3)
  #  if (LC=='LC1') {F_eff.LC[,,h] = (A_r*f_eff.r^3 + A_veg*SSR[,,h]^3 + A_mix*(SSR[,,h]*f_eff.r)^3)^(1/3)}
  #  else if (LC=='LC0') {F_eff.LC[,,h] = (frc.r.C1*f_eff.r^3 + (frc.veg.C1+frc.mix.C1)*SSR[,,h]^3)^(1/3)}  # newly define LC0, removed mixed regime
  #}
  if (LC=='LC1') {
    f_eff_veg.comp = sweep(SSR^3, c(1,2), (frc.veg.C1), '*')
    f_eff_mix.comp = sweep(SSR^3, c(1,2), (frc.mix.C1*(f_eff.r^3)), '*')
    F_eff.LC = (sweep(f_eff_veg.comp+f_eff_mix.comp, c(1,2), (frc.r.M2*f_eff.r^3), '+')) ^ (1/3)
    rm('f_eff_veg.comp','f_eff_mix.comp')  
  } else if (LC=='LC0') {
    f_eff_veg.comp = sweep(SSR^3, c(1,2), (frc.veg.C1+frc.mix.C1), '*')
    F_eff.LC = (sweep(f_eff_veg.comp, c(1,2), (frc.r.C1*f_eff.r^3), '+')) ^ (1/3)
    rm(f_eff_veg.comp)
  }
  
  
  #plot.field(F_eff.LC[,,6], LON.C1, LAT.C1, type='def', zlim=c(0,1), col=WBGYR)
  
  # then get hybrid aeolian roughness length (Z0.LC1) from hybrid F_eff, using Marticorena's equation
  #f_eff.r = 1 - ( log(z0.r/z0s) / log(0.7*(X/z0s)^0.8) )  
  #Z0.LC = z0s * exp( (1 - F_eff.LC) * log(0.7*(X/z0s)^0.8) )   # in m
  #plot.field.log(Z0.LC1[,,6], LON.C1, LAT.C1, type='def', zlim=c(-5,0), col=tim.colors(32)[32:1], legend.mar=5)
  #plot.field.log(Z0M[,,6], LON.C1, LAT.C1, type='def', zlim=c(-5,0), col=tim.colors(32)[32:1], legend.mar=5)
  #plot.field.log(z0.r*landfilt.C1, LON.C1, LAT.C1, type='def', zlim=c(-5,0), col=tim.colors(32)[32:1], legend.mar=5)   # also in m
  
  #print('finished drag partition calculation; now select friction velocity at relevant scale')
  
  ###################
  # reconstruct friction velocity with different meanings and choose which to use
  
  # test if U10, USTAR, and Z0M are consistent with each other
  
  #for (aa in 1) {
  if (wind.scale=='default') {print('default MERRA2 ustar is used')}
  else {
    Z0.LC = z0s * exp( (1 - F_eff.LC) * log(0.7*(X/z0s)^0.8) )   # in m
    U10   = .env$U10.C1[i1:il,j1:jl,]  # m/s. Define U10 as the total wind speed of U10M and V10M. Don't confuse U10 with U10M. 
    if (wind.scale=='aeo_nonneu') {
      KSI_m = k*U10/USTAR - log(10 / Z0M) # first determine MOST instability (KSI)
      USTAR.aeo.nonneu = k * U10 / (log(10 / Z0.LC)+KSI_m); print('aeolian+non-neutral ustar is used')}  # reconstruct USTAR from U10 and Z0M
    #USTAR.aero.neu = k * U10 / log(10 / Z0M)  # reconstruct USTAR from U10 and Z0M; aero means aerodynamic, neu means neutral (log wind profile) 
    # calculate aeolian ustar assuming neutral condition
    else if (wind.scale=='aeo_neu') {USTAR.aeo.neu = k * U10 / log(10 / Z0.LC); print('aeolian+neutral ustar is used')}  # reconstruct USTAR from U10 and Z0M
  }
  #plot.field(USTAR.aero.neu[,,8] / USTAR[,,8], LON.C1, LAT.C1, type='def', zlim=c(0,2), col=TEMP_DIFF_65)
  # then also reconstruct ustar under non-neutral condition
  #plot.field(USTAR[,,6]*landfilt.C1, LON.C1, LAT.C1, type='def', zlim=c(0,0.8), col=WBGYR)
  #plot.field(USTAR.aeo.nonneu[,,6], LON.C1, LAT.C1, type='def', zlim=c(0,0.8), col=WBGYR)
  #plot.field(USTAR.aeo.nonneu[,,6]/USTAR.aeo.neu[,,6], LON.C1, LAT.C1, type='def', zlim=c(0,2), col=TEMP_DIFF_65)
  #plot.field(USTAR.aeo.neu[,,6], LON.C1, LAT.C1, type='def', zlim=c(0,0.8), col=WBGYR)
  #A = apply(USTAR, c(1,2), mean, na.rm=T)
  #plot.field(A*landfilt.C1, LON.C1, LAT.C1, legend.mar=5, type='def', zlim=c(0,0.7))
  #A = apply(USTAR.aeo.nonneu, c(1,2), mean, na.rm=T)
  #plot.field(A, LON.C1, LAT.C1, legend.mar=5, type='def', zlim=c(0,0.7))
  #A = apply(USTAR.aeo.neu, c(1,2), mean, na.rm=T)
  #plot.field(A, LON.C1, LAT.C1, legend.mar=5, type='def', zlim=c(0,0.7))
  
  #print('finished threshold calculation; now calculate fluid and impact thresholds')
  
  ###################
  # construct dust emission thresholds
  
  # first get dry fluid threshold
  ustar_ft.dry = 1.0*sqrt(0.0123 * (dns_slt*g*D_p + gamma/D_p)) / sqrt(RHOA)
  #plot.field(ustar_ft.dry[,,6], LON.C1, LAT.C1)
  
  # to get wet fluid threshold,
  # get gravimetric water content from volumetric water content
  gwc_sfc = SHR_CONST_RHOFW*sweep(SFMC, c(1,2), bulk_den, '/')
  #plot.field(gwc_sfc[,,6], LON.C1, LAT.C1)
  #plot.field.log(gwc_sfc[,,6], LON.C1, LAT.C1, type='def',zlim=c(-3,0), legend.mar=5)
  
  # then get soil moisture effect on fluid threshold
  #ptm = proc.time()
  #frc_thr_wet_fct = array(NaN, dim=c(length(indi), length(indj), length(indt)))
  #for (h in 1:length(indt)) {
  #  A = sqrt(1.0 + 1.21*(100.0*(gwc_sfc[,,h] - gwc_thr))^0.68)#  A[which(gwc_sfc[,,h] < gwc_thr)] = 1
  #  frc_thr_wet_fct[,,h] = A
  #}
  #proc.time() - ptm
  #rm(A)  # remove dummy
  #ptm = proc.time()
  frc_thr_wet_fct = sqrt(1.0 + 1.21*(100.0* (sweep(gwc_sfc, c(1,2), gwc_thr, '-')) )^0.68)
  frc_thr_wet_fct[which(sweep(gwc_sfc, c(1,2), gwc_thr, '<'))] = 1
  #proc.time() - ptm
  #plot.field(frc_thr_wet_fct[,,6], LON.C1, LAT.C1, type='def', zlim=c(1,3.6))
  
  # then get fluid and impact thresholds (now assume using FAO's moisture effect)
  ustar_ft.wet = ustar_ft.dry * frc_thr_wet_fct
  #ustar_it = ustar_ft.dry * B_it * landfilt.C1
  ustar_it = B_it * sweep((ustar_ft.dry), c(1,2), landfilt.C1, '*')
  #plot.field(ustar_it[,,6], LON.C1, LAT.C1)
  #plot.field(ustar_ft.wet[,,6], LON.C1, LAT.C1)
  #A = apply(USTAR - ustar_it, c(1,2), mean, na.rm=T)
  #plot.field(A*landfilt.C1, LON.C1, LAT.C1, legend.mar=5, type='sign', zlim=c(-0.5,0.5), col=TEMP_DIFF_65)
  #A = apply(USTAR.aeo.nonneu - ustar_it, c(1,2), mean, na.rm=T)
  #plot.field(A, LON.C1, LAT.C1, legend.mar=5, type='sign', zlim=c(-0.5,0.5), col=TEMP_DIFF_65)
  
  #print('finished threshold calculation; now calculate intermittency effect')
  
  ###################
  # construct intermittency
  
  
  # choose a friction velocity. Aeolian/Aerodynamic? Neutral/Non-neutral?
  # USTAR.aeo.neu, USTAR.aeo.nonneu, USTAR.aero.neu, USTAR (aero+nonneu, default MERRA2 output)
  for (aa in 1) {
    if (wind.scale=='default') {wnd_frc_slt = USTAR*F_eff.LC; print('default MERRA2 ustar is used for intermittency and emission')}   #  default MERRA2 (aero+nonneu)output
    else if (wind.scale=='aeo+neu') {wnd_frc_slt = USTAR.aeo.neu*F_eff.LC; print('aeolian+neutral ustar is used for intermittency and emission')}   #  reconstructed aeo+neu USTAR
    else if (wind.scale=='aeo+nonneu') {wnd_frc_slt = USTAR.aeo.nonneu*F_eff.LC; print('aeolian+neutral ustar is used for intermittency and emission')}   #  reconstructed aeo+nonneu USTAR
  }
  
  if (use.intermittency=='TRUE') {
    # get Monin-Obukhov length in m (Donald Golder, 1972; Khaled Essa, 1999)
    # we use aerodynamic non-neutral USTAR here
    Obu_L = - RHOA * c_p * T10M * USTAR^3 / (k * g * SHLAND)
    
    u_mean_slt = (wnd_frc_slt/k) * log(0.1 / 1e-4)  # We used z0a = 1e-4 m in C1 and assumed. Should we change to aeolian data here?
    # sd of lowpass-filtered wind speed
    # u_sd_slt = wnd_frc_slt * (12 )^0.333   # note that there is no buoyancy term here; need to add in later when Obukhov length is obtained (1 Oct 2021)
    u_sd_slt = wnd_frc_slt * (12 - 0.5 * PBLH/Obu_L)^0.333   # buoyancy term is added
    # then get thresholds at saltation height
    u_fld_thr = (ustar_ft.wet/k) * log(0.1 / 1e-4) # to avoid model error
    u_impct_thr = (ustar_it/k) * log(0.1 / 1e-4)  # to avoid model error
    # threshold crossing rate
    thr_crs_rate = (exp((u_fld_thr^2 - u_impct_thr^2 - 2 * u_mean_slt * (u_fld_thr - u_impct_thr)) / (2 * u_sd_slt^2)) + 1)^(-1)
    # probability that lowpass-filtered wind speed does not exceed u_ft
    prb_crs_fld_thr = 0.5 * (1 + erf((u_fld_thr - u_mean_slt) / (1.414 * u_sd_slt)))
    # probability that lowpass-filtered wind speed does not exceed u_it
    prb_crs_impct_thr = 0.5 * (1 + erf((u_impct_thr - u_mean_slt) / (1.414 * u_sd_slt)))
    # intermittency factor (from 0 to 1)
    intrmtncy_fct = 1 - prb_crs_fld_thr + thr_crs_rate * (prb_crs_fld_thr - prb_crs_impct_thr)
  }
  #plot.field(u_mean_slt[,,6], LON.C1, LAT.C1)
  #plot.field(u_sd_slt[,,6], LON.C1, LAT.C1)
  #plot.field(u_fld_thr[,,6], LON.C1, LAT.C1)
  #plot.field(u_impct_thr[,,6], LON.C1, LAT.C1)
  #plot.field(intrmtncy_fct[,,6], LON.C1, LAT.C1)
  #A = apply(intrmtncy_fct, c(1,2), mean, na.rm=T)
  #plot.field(A, LON.C1, LAT.C1, col=WBGYR, legend.mar=5, type='def', zlim=c(0,1))
  
  #print('finished intermittency; now calculate dust emission flux')
  
  ###################
  # a few more essential variables for dust emission equations
  # standardized fluid threshold (m/s)
  ustar_st = ustar_ft.wet * sqrt(RHOA/rho_a0)
  # dust emission coefficient (a measure of soil erodibility)
  C_d = C_d0 * exp(-C_e * (ustar_st-ustar_st0)/ustar_st0);
  #A = apply(C_d, c(1,2), mean, na.rm=T)
  #plot.field.log(A/C_d0, LON.C1, LAT.C1, type='def', zlim=c(-5,0), legend.mar=5)
  
  
  ###################
  # compute 3-D dust emissions
  
  # first get fragmentation exponent and limit to frag_expt_lim and below
  frag_expt = (C_alpha*(ustar_st-ustar_st0)/ustar_st0)
  frag_expt[which(frag_expt>frag_expt_lim)] = frag_expt_lim
  #A = apply(frag_expt, c(1,2), mean, na.rm=T)
  #plot.field(A, LON.C1, LAT.C1, legend.mar=5)
  
  # compute dust emissions
  for (aa in 1) {
    if (use.intermittency == 'TRUE') {
      F_d = array(0, dim=dim(C_d))   # using dimension of C_d to define F_d dimension
      ind = which(wnd_frc_slt > ustar_it)  # 3-D index for which wind speed > threshold
      F_d[ind] = C_tune * C_d[ind] * f_bare[ind] * RHOA[ind] * (wnd_frc_slt[ind]^2-ustar_it[ind]^2)/ustar_it[ind] * (wnd_frc_slt[ind]/ustar_it[ind])^(frag_expt[ind]) * intrmtncy_fct[ind]             # calculate F_d for indices with wind speed > threshold
      F_d = sweep(F_d, c(1,2), mss_frc_cly_vld, '*')  # multiply 3-D fields by 2-D fields
      #F_d_sum = apply(F_d, c(1,2), sum, na.rm=T) * 3600 * dt * 365
      #plot.field.log(F_d_sum, LON.C1, LAT.C1, legend.mar=5, col=WBGYR)
    } else if (use.intermittency == 'FALSE') {
      F_d = array(0, dim=dim(C_d))   # using dimension of C_d to define F_d dimension
      ind = which(wnd_frc_slt > ustar_ft.wet)  # 3-D index for which wind speed > threshold
      F_d[ind] = C_tune * C_d[ind] * f_bare[ind] * RHOA[ind] * (wnd_frc_slt[ind]^2-ustar_ft.wet[ind]^2)/ustar_ft.wet[ind] * (wnd_frc_slt[ind]/ustar_ft.wet[ind])^(frag_expt[ind])              # calculate F_d for indices with wind speed > threshold
      F_d = sweep(F_d, c(1,2), mss_frc_cly_vld, '*')  # multiply 3-D fields by 2-D fields
    }
  }
  
  #filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=3/dust_emis_flx_1x1_', date_vec[d], '.RData', sep='')
  filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/dust_emis_flx_1x1_', date_vec[d], '.RData', sep='')
  save('F_d', 'LON.C1', 'LAT.C1', file=filename)
  # 'C_d','f_bare','frag_expt','KSI_m' don't change with roughness and drag partition factor, so we save them for easier computations in upscaling process.
  #filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=3/vars_for_upscale_1x1_', date_vec[d], '.RData', sep='')
  #save('C_d','f_bare','frag_expt','KSI_m', file=filename)
  
  print('emission and other variables saved')
  
}




ustar_it.save = ustar_ft.wet.save = array(0, dim=c(length(LON.C1), length(LAT.C1), 12))
for (d in 1:length(date_vec)) {# d-th day
  print(paste('day ', d, ': ', date_vec[d], sep=''))
  #print('extract MERRA2 data', sep='')
  ###################
  # get 3D data in time and space, subsetted as regional data with variable time resolutions
  i1 = indi[1]; j1 = indj[1]; t1 = indt[1]
  il = length(indi); jl = j1+length(indj)-1; tl = t1+length(indt)
  
  #ptm = proc.time()
  #nc = nc_open(files_flx[d])#; names(nc$var)#, start=c(11,26), count=c(5,5)
  load(files[d], envir=.env)#; ls(.env)
  RHOA  = .env$RHOA.C1[i1:il,j1:jl,]     # kg/m3
  #proc.time() - ptm
  SFMC  = .env$SFMC.C1[i1:il,j1:jl,]
  
  #print('finished loading MERRA2 met fields; now calculate hybrid drag partition and hybrid roughness length')
  
  ###################
  # construct dust emission thresholds
  
  # first get dry fluid threshold
  ustar_ft.dry = 1.0*sqrt(0.0123 * (dns_slt*g*D_p + gamma/D_p)) / sqrt(RHOA)
  #plot.field(ustar_ft.dry[,,6], LON.C1, LAT.C1)
  
  # to get wet fluid threshold,
  # get gravimetric water content from volumetric water content
  gwc_sfc = SHR_CONST_RHOFW*sweep(SFMC, c(1,2), bulk_den, '/')
  
  # then get soil moisture effect on fluid threshold
  frc_thr_wet_fct = sqrt(1.0 + 1.21*(100.0* (sweep(gwc_sfc, c(1,2), gwc_thr, '-')) )^0.68)
  frc_thr_wet_fct[which(sweep(gwc_sfc, c(1,2), gwc_thr, '<'))] = 1
  
  # then get fluid and impact thresholds (now assume using FAO's moisture effect)
  ustar_ft.wet = ustar_ft.dry * frc_thr_wet_fct
  ustar_ft.wet.save = ustar_ft.wet.save + ustar_ft.wet
  #ustar_it = ustar_ft.dry * B_it * landfilt.C1
  ustar_it = B_it * sweep((ustar_ft.dry), c(1,2), landfilt.C1, '*')
  ustar_it.save = ustar_it.save + ustar_it
  
}





ustar_ft.wet.save.yr = apply(ustar_ft.wet.save, c(1,2), mean, na.rm=T)/length(date_vec)
ustar_it.save.yr = apply(ustar_it.save, c(1,2), mean, na.rm=T)/length(date_vec)

plot.field.proj(ustar_ft.wet.save.yr, LON.C1, LAT.C1, col=WBGYR)
plot.field(ustar_ft.wet.save.yr, LON.C1, LAT.C1, col=WBGYR, legend.only=T, zlim=c(0.27,1.4), legend.width=0.8, horizontal=T)
plot.field.proj(ustar_it.save.yr, LON.C1, LAT.C1, col=WBGYR)


ustar_ft.wet.save.yr1 = apply(ustar_ft.wet.save, c(1,2), mean, na.rm=T)/length(date_vec)
plot.field.proj(ustar_ft.wet.save.yr1/ustar_ft.wet.save.yr, LON.C1, LAT.C1, col=WBGYR)
