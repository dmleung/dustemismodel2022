#############################################################
#############################################################

# dust emission code version 1: 2006 MERRA2 0.5x0.625 data for dust emissions
# extract N Africa for calculation only (or specify other spatial domains)
# allow 1-hourly, 2-hourly data or other time resolution simulations (can skip some timestep and multiply dt at the end if one think e.g. 2-hourly sampling is good enough)
# allow showing actual dust emissions and annually scaled emissions (Tg / yr)
# set up spatial and temporal domain

# 30 Sep 2021
# remaining issues:
# 1. Confirm does CESM using porosity or effective porosity (discounting ice portion) to calculate gravimetric water content? We use MERRA2 here which uses porosity
# 2. CESM uses liquid water fractions for dust emissions, but here we haven't added this in for MERRA2 dust emissions
# 3. For intermittency, we assumed z0a = 1e-4 m globally. Change to hybrid Z0a? Are we going to do the same for CESM?
# MERRA2 processing trial
# 25 Sep 2021
# Main script for the first project
# 6 Oct 2021
# Now I a big time loop for everything; but if later on find that the loop takes too long, one can easily change the script to several smaller time loops and run one by one.
# 12 Oct 2021
# Horn of Africa has the largest emissions across the globe, and also NE China has the second largest emissions. Any ways to suppress those emissions (e.g., further reduce the fragmentation exponent?)
# 18 Nov 2021
# The code has been revised to include Zender et al. (2003) as the dust production model (DPM). Other codes mostly only have Kok et al. (2014) as the DPM.

# set working directory
setwd('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR')
#setwd("/Volumes/GoogleDrive-103510611982741156765/My Drive/CESM/CESMDUSTinR")
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
#indi = 251:390; indj = 171:260     # lon and lat indices for N Africa
indi = 1:576; indj = 65:350 
# choose clay dataset: 'FAO' or 'SG' (SoilGrids)
clay_dataset = 'FAO'
# choose the wind scale for dust emission, including default (aerodynamic + nonneutral) ustar (default), aeolian + nonneutral star (aeo+nonneu), and aeolian + neutral ustar (aeo+neu), and default + no drag partition (no.drag.partition)
wind.scale='default'
# choose land cover averaging method: LC0 (rock and veg) or LC1 (rock, veg and mix)
LC = 'LC0'
# choose whether default Kok et al. (2014) scheme (K14) or intermittency (K14+Comola et al., 2019) is used. 'FALSE' is K14, 'TRUE' is K14+C19 intermittency scheme.
use.intermittency = TRUE
# if use.intermittency is FALSE, need to determine whether one wants to fluid or impact threshold; if use.intermittency is TRUE then this parameter is not important
use.threshold = 'impact'   
# fragmentation exponent limit (CESM simulations can blow up when frag_expt is too large when using impact threshold for dust emission equation)
frag_expt_lim = 3   
# dust production model (DPM): Kok et al. (2014) K14 or Zender et al. (2003) Z03, added 18 Nov 2021
DPM = 'K14'
if (DPM=='Z03') use.intermittency = FALSE # by default Z03 has no intermittency eff
# if source function is needed (for Z03), choose 'Zender' or 'Ginoux'. For K14, type 'FALSE'
src_fn = 'FALSE'
#directory name ('K14_default', 'K14_u_ft_Dp=130', 'K14_u_ft_Dp=130_F_eff', 'K14_u_it_Dp=130', 'K14_u_it_Dp=130_F_eff', 'L21_final')
directory = 'L21_final'
##################
# calculate model parameters based on user inputs
date_vec = make.date.vec(startdate, enddate)  # make date vec
time = 1:(length(date_vec)*24/dt)  # a vector of timesteps
indt = seq(dt, 24, by=dt)   # hour of day
# get lon lat
nc = nc_open('../../MERRA2/upscale_2006/MERRA2_300.tavg1_2d_lnd_Nx.20060101.SUB.nc')
lon.M2 = ncvar_get(nc, 'lon')[indi]; lat.M2 = ncvar_get(nc, 'lat')[indj]
nc_close(nc)

###################
# constant parameters
C_tune = 0.05   # global tuning constant
gamma = 1.65e-4;   # soil aggregation effect on dust emission threshold
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
cst_slt = 2.61 # saltation constant for Z03, added 18 Nov 2021
flx_mss_fdg_fct = 5.0e-4 # Empirical emission tuning factor for Z03, added 18 Nov 2021
###################
# load required time-invariant datasets
.env = new.env()   # set a new environment for storing data temporarily
# load GLCNMO LULC data
load('GLCNMO_frc_area_05x0625.RData',envir=.env); ls(.env)
frc.r.M2 = .env$frc.r.M2[indi,indj]; frc.veg.M2 = .env$frc.veg.M2[indi,indj]
frc.mix.M2 = .env$frc.mix.M2[indi,indj]
#A_r = frc.r.M2 / (frc.r.M2+frc.veg.M2+frc.mix.M2)
#A_veg = frc.veg.M2 / (frc.r.M2+frc.veg.M2+frc.mix.M2)
#A_mix = frc.mix.M2 / (frc.r.M2+frc.veg.M2+frc.mix.M2)
#rm('frc.r.M2','frc.veg.M2','frc.mix.M2')
# load landfilt 
load('landfilt_05x0625.RData',envir=.env)
landfilt.M2 = .env$landfilt.M2[indi,indj]
# get Prigent 2005 roughness length and rock drag partition factor
load('Pr05_z0.min_05x0625.RData',envir=.env)
z0.r = .env$z0.M2[indi,indj]/100   # in cm, divide by 100 to convert to m
# in this formulation, set globally constant z0s because soil median diameter ~ 130 um
f_eff.r = 1 - ( log(z0.r/z0s) / log(0.7*(X/z0s)^0.8) )  
f_eff.r = f_eff.r*landfilt.M2   # take away values in Caspian sea
# image(f_eff.r)
# get porosity  (issue: MERRA2 porosity may be different from CESM's effective porosity)
nc = nc_open('/Volumes/GoogleDrive/My Drive/MERRA2/upscale_2006/MERRA2_100.const_2d_lnd_Nx.00000000.nc4')#; names(nc$var)
poros = ncvar_get(nc, 'poros')[indi,indj] # porosity (m3 / m3, or volumetric soil water) 
#image(poros)
#plot.field(poros, lon.M2, lat.M2)
# get bulk density of dry soil material [kg/m^3]
bulk_den = (1 - poros)*dns_slt 
#plot.field(bulk_den, lon.M2, lat.M2, legend.mar=4)
# use CESM's FAO clay fraction
for (aa in 1) {
  if (clay_dataset == 'FAO') {
    load('f_clay_FAO_05x06.RData',envir=.env)
    f_clay.FAO.05x06 = .env$f_clay.FAO.05x06[indi,indj]
    #plot.field(f_clay.FAO.05x06, lon.M2, lat.M2, type='def', zlim=c(0,0.65))
    # threshold gravimetric moisture and clay limitation on dust emission
    gwc_thr = 0.01*(0.17*f_clay.FAO.05x06*100 + 0.0014*(f_clay.FAO.05x06*100)^2)
    #plot.field(gwc_thr, lon.M2, lat.M2)
    mss_frc_cly_vld = f_clay.FAO.05x06
    mss_frc_cly_vld[which(mss_frc_cly_vld>0.2)] = 0.2
    #plot.field(mss_frc_cly_vld, lon.M2, lat.M2)
    print('used FAO clay data')
  }  else if (clay_dataset == 'SG') {
    load('f_clay_SoilGrids_05x06.RData',envir=.env)
    f_clay.SG.05x06 = .env$f_clay.SG.05x06[indi,indj]
    #plot.field(f_clay.SG.05x06, lon.M2, lat.M2, type='def', zlim=c(0,0.65))
    # threshold gravimetric moisture and clay limitation on dust emission
    gwc_thr = 0.01*(0.17*f_clay.SG.05x06*100 + 0.0014*(f_clay.SG.05x06*100)^2)
    #plot.field(gwc_thr, lon.M2, lat.M2)
    mss_frc_cly_vld = f_clay.SG.05x06
    mss_frc_cly_vld[which(mss_frc_cly_vld>0.2)] = 0.2
    #plot.field(mss_frc_cly_vld, lon.M2, lat.M2)
    print('used SoilGrids clay data')
  }
}
# calculate sandblasting efficiency if Z03 is used, added 18 Nov 2021
for (aa in 1) {
  if (DPM=='Z03') {
    # sandblasting efficiency
    dst_slt_flx_rat_ttl = 100*exp( log(10) * (13.4*mss_frc_cly_vld - 6))
  }
  #plot.field(dst_slt_flx_rat_ttl, lon.M2, lat.M2)
}
# choose source function if using Z03 dust emission equation, added 18 Nov 2021
for (aa in 1) {
  if (src_fn=='Ginoux') {
    load('ginoux_dust_source_05x0625.RData',envir=.env)
    mbl_bsn_fct = .env$mbl_bsn_fct_G.05x06[indi,indj]
    print('used Ginoux source function')
  } else if (src_fn=='Zender') {
    load('zender_dust_source_05x0625.RData',envir=.env)
    mbl_bsn_fct = .env$mbl_bsn_fct_Z.05x06[indi,indj]
    print('used Zender source function')
  }
}
#plot.field.log(mbl_bsn_fct, lon.M2, lat.M2, legend.mar=5)
###################
# scrap files for selected dates
files_flx = files_slv = files_lnd = NULL   # character strings storing the file names
for (dd in 1:length(date_vec)) {
  date_vec[dd]
  #files_flx = c(files_flx, Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/tavg1_2d_flx_Nx/MERRA2*flx*', date_vec[dd], '*.nc', sep='')))  # for surface flux variables: surface USTAR, Z0M, and RHOA (air density)
  #files_slv = c(files_slv, Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/tavg1_2d_slv_Nx/MERRA2*slv*', date_vec[dd], '*.nc', sep='')))  # for surface level variables: surface U10, V10, and DISPH (displacement height)
  #files_lnd = c(files_lnd, Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/tavg1_2d_lnd_Nx/MERRA2*lnd*', date_vec[dd], '*.nc', sep='')))  # for land variables: SFMC (Volumetric water content) and LAI
  files_flx = c(files_flx, Sys.glob(paste('/Volumes/GoogleDrive/My Drive/upscale_2006/tavg1_2d_flx_Nx/MERRA2*flx*', date_vec[dd], '*.nc', sep='')))  # for surface flux variables: surface USTAR, Z0M, and RHOA (air density)
  files_slv = c(files_slv, Sys.glob(paste('/Volumes/GoogleDrive/My Drive/upscale_2006/tavg1_2d_slv_Nx/MERRA2*slv*', date_vec[dd], '*.nc', sep='')))  # for surface level variables: surface U10, V10, and DISPH (displacement height)
  files_lnd = c(files_lnd, Sys.glob(paste('/Volumes/GoogleDrive/My Drive/upscale_2006/tavg1_2d_lnd_Nx/MERRA2*lnd*', date_vec[dd], '*.nc', sep='')))  # for land variables: SFMC (Volumetric water content) and LAI
}

###################
# start the main loop
# for days in d
# d = 1   # d-th day

ptm = proc.time()
#for (d in 1) {# d-th day
for (d in 1:length(date_vec)) {# d-th day
  
print(paste('day ', d, ': ', date_vec[d], sep=''))
#print('extract MERRA2 data', sep='')
###################
# get 3D data in time and space, subsetted as regional data with variable time resolutions
i1 = indi[1]; j1 = indj[1]; #t1 = indt[1]
il = length(indi); jl = length(indj); #tl = length(indt)

#ptm = proc.time()
#nc = nc_open(files_flx[d]); names(nc$var)#, start=c(11,26), count=c(5,5)
#USTAR = ncvar_get(nc, 'USTAR')[indi,indj,indt]   # m/s
#Z0M = ncvar_get(nc, 'Z0M')[indi,indj,indt]       # m
#RHOA = ncvar_get(nc, 'RHOA')[indi,indj,indt]     # kg/m3
#PBLH = ncvar_get(nc, 'PBLH')[indi,indj,indt]     # kg/m3
#nc = nc_open(files_slv[d]); names(nc$var)
#DISPH = ncvar_get(nc, 'DISPH')[indi,indj,indt]   # m
#U10M = ncvar_get(nc, 'U10M')[indi,indj,indt]     # m/s
#V10M = ncvar_get(nc, 'V10M')[indi,indj,indt]     # m/s
#T10M = ncvar_get(nc, 'T10M')[indi,indj,indt]     # K
#U10 = sqrt(U10M^2 + V10M^2)  # m/s. Define U10 as the total wind speed of U10M and V10M. Don't confuse U10 with U10M. 
#proc.time() - ptm

#ptm = proc.time()
nc = nc_open(files_flx[d])#; names(nc$var)#, start=c(11,26), count=c(5,5)
USTAR = ncvar_get(nc, 'USTAR', start=c(i1,j1,1), count=c(il,jl,24)); USTAR = USTAR[,,indt] # m/s
Z0M = ncvar_get(nc, 'Z0M', start=c(i1,j1,1), count=c(il,jl,24)); Z0M = Z0M[,,indt]     # m
RHOA = ncvar_get(nc, 'RHOA', start=c(i1,j1,1), count=c(il,jl,24)); RHOA = RHOA[,,indt]     # kg/m3
PBLH = ncvar_get(nc, 'PBLH', start=c(i1,j1,1), count=c(il,jl,24)); PBLH = PBLH[,,indt]     # kg/m3
nc = nc_open(files_slv[d])#; names(nc$var)
#DISPH = ncvar_get(nc, 'DISPH', start=c(i1,j1,1), count=c(il,jl,24))[,,indt]   # m
#U10M = ncvar_get(nc, 'U10M', start=c(i1,j1,1), count=c(il,jl,24)); U10M = U10M[,,indt]     # m/s
#V10M = ncvar_get(nc, 'V10M', start=c(i1,j1,1), count=c(il,jl,24)); V10M = V10M[,,indt]     # m/s
T10M = ncvar_get(nc, 'T10M', start=c(i1,j1,1), count=c(il,jl,24)); T10M = T10M[,,indt]     # K
#U10 = sqrt(U10M^2 + V10M^2)  # m/s. Define U10 as the total wind speed of U10M and V10M. Don't confuse U10 with U10M. 
#proc.time() - ptm
#rm('U10M', 'V10M') # remove data to save space
nc = nc_open(files_lnd[d])#; names(nc$var)
SFMC = ncvar_get(nc, 'SFMC', start=c(i1,j1,1), count=c(il,jl,24)); SFMC = SFMC[,,indt]
LAI = ncvar_get(nc, 'LAI', start=c(i1,j1,1), count=c(il,jl,24)); LAI = LAI[,,indt]
SHLAND = ncvar_get(nc, 'SHLAND', start=c(i1,j1,1), count=c(il,jl,24)); SHLAND = SHLAND[,,indt]
#plot.field(DISPH[,,3], lon.M2, lat.M2)
#plot.field(USTAR[,,3], lon.M2, lat.M2)
#plot.field(Z0M[,,3], lon.M2, lat.M2)
#plot.map.log((Z0M[,,3]),lon.M2, lat.M2,type='def',zlim=c(-4,0), legend.mar=5)

#print('finished loading MERRA2 met fields; now calculate hybrid drag partition and hybrid roughness length')
#proc.time() - ptm
###################
# get drag partition factors and hybrid roughness length (later on, make this as a separate time loop independent of dust emission calculation?)

# get vegetation roughness length following Okin's scheme using LAI
# first use LAI to get bare land fraction using Natalie Mahowald's equation
# and vegetation drag partition using Okin-Pierre parameterization
# Natalie's bare land fraction approximation
f_bare = 1 - LAI/LAI_thr; f_bare[which(f_bare<0)] = NaN
#f_bare_0.3 = 1 - LAI[,,3]/0.3; f_bare_0.3[which(f_bare_0.3<0)] = NaN
#plot.field(f_bare[,,6], lon.M2, lat.M2)
#plot.field(f_bare_0.3, lon.M2, lat.M2)
#plot.field.log(LAI[,,6], lon.M2, lat.M2, zlim=c(-2,1), legend.mar=5)

# make a filter for dust emission sources
#emisfilt = f_bare/f_bare    # this filter is 3D and changes with time; or one can make a annual, 2D one. 
#image(emisfilt)

# Okin's vegetation drag partition
# Using Caroline Pierre (2014)'s formulation of Okin's scheme
K = pi/2 * (1/LAI - 1); K[which(K<0)] = 0
SSR = (K + U_0 * c) / (K + c)  # Pierre-Okin shear stress ratio (SSR)
#plot.field(SSR[,,3], lon.M2, lat.M2, type='def', zlim=c(0.6,0.9))

# combine rock and vegetation drag partition factors
# then combine drag partition factors
#F_eff.LC = array(NaN, dim=c(length(indi), length(indj), length(indt)))
#F_eff.LC2 = array(NaN, dim=c(length(indi), length(indj), length(indt)))
#for (h in 1:length(indt)) {
  #F_eff.LC1[,,h] = (frc.r.M2*f_eff.r^3 + frc.veg.M2*SSR[,,h]^3 + frc.mix.M2*(SSR[,,h]*f_eff.r)^3)^(1/3)
#  if (LC=='LC1') {F_eff.LC[,,h] = (frc.r.M2*f_eff.r^3 + frc.veg.M2*SSR[,,h]^3 + frc.mix.M2*(SSR[,,h]*f_eff.r)^3)^(1/3)}
#  else if (LC=='LC0') {F_eff.LC[,,h] = (frc.r.M2*f_eff.r^3 + (frc.veg.M2+frc.mix.M2)*SSR[,,h]^3)^(1/3)}  # newly define LC0, removed mixed regime
#}
if (LC=='LC1') {
  f_eff_veg.comp = sweep(SSR^3, c(1,2), (frc.veg.M2), '*')
  f_eff_mix.comp = sweep(SSR^3, c(1,2), (frc.mix.M2*(f_eff.r^3)), '*')
  F_eff.LC = (sweep(f_eff_veg.comp+f_eff_mix.comp, c(1,2), (frc.r.M2*f_eff.r^3), '+')) ^ (1/3)
  rm('f_eff_veg.comp','f_eff_mix.comp')  
} else if (LC=='LC0') {
  f_eff_veg.comp = sweep(SSR^3, c(1,2), (frc.veg.M2+frc.mix.M2), '*')
  F_eff.LC = (sweep(f_eff_veg.comp, c(1,2), (frc.r.M2*f_eff.r^3), '+')) ^ (1/3)
  rm(f_eff_veg.comp)
}

#plot.field(F_eff.LC1[,,6], lon.M2, lat.M2, type='def', zlim=c(0,1), col=WBGYR)
#plot.field(F_eff.LC2[,,6], lon.M2, lat.M2, type='def', zlim=c(0,1), col=WBGYR)

# then get hybrid aeolian roughness length (Z0.LC1) from hybrid F_eff, using Marticorena's equation
#f_eff.r = 1 - ( log(z0.r/z0s) / log(0.7*(X/z0s)^0.8) )  
#Z0.LC = z0s * exp( (1 - F_eff.LC) * log(0.7*(X/z0s)^0.8) )   # in m
#plot.field.log(Z0.LC1[,,6], lon.M2, lat.M2, type='def', zlim=c(-5,0), col=tim.colors(32)[32:1], legend.mar=5)
#plot.field.log(Z0M[,,6], lon.M2, lat.M2, type='def', zlim=c(-5,0), col=tim.colors(32)[32:1], legend.mar=5)
#plot.field.log(z0.r*landfilt.M2, lon.M2, lat.M2, type='def', zlim=c(-5,0), col=tim.colors(32)[32:1], legend.mar=5)   # also in m

#print('finished drag partition calculation; now select friction velocity at relevant scale')

###################
# reconstruct friction velocity with different meanings and choose which to use

# test if U10, USTAR, and Z0M are consistent with each other

#for (aa in 1) {
if (wind.scale=='default') {print('default MERRA2 ustar is used')}
else {
  Z0.LC = z0s * exp( (1 - F_eff.LC) * log(0.7*(X/z0s)^0.8) )   # in m
  nc = nc_open(files_slv[d])#; names(nc$var)
  U10M = ncvar_get(nc, 'U10M', start=c(i1,j1,1), count=c(il,jl,24)); U10M = U10M[,,indt]     # m/s
  V10M = ncvar_get(nc, 'V10M', start=c(i1,j1,1), count=c(il,jl,24)); V10M = V10M[,,indt]     # m/s
  U10 = sqrt(U10M^2 + V10M^2)  # m/s. Define U10 as the total wind speed of U10M and V10M. Don't confuse U10 with U10M. 
  rm('U10M', 'V10M') # remove data to save space
    
  if (wind.scale=='aeo_nonneu') {
    KSI_m = k*U10/USTAR - log(10 / Z0M) # first determine MOST instability (KSI)
    USTAR.aeo.nonneu = k * U10 / (log(10 / Z0.LC)+KSI_m); print('aeolian+non-neutral ustar is used')}  # reconstruct USTAR from U10 and Z0M
  #USTAR.aero.neu = k * U10 / log(10 / Z0M)  # reconstruct USTAR from U10 and Z0M; aero means aerodynamic, neu means neutral (log wind profile) 
  # calculate aeolian ustar assuming neutral condition
  else if (wind.scale=='aeo_neu') {
    USTAR.aeo.neu = k * U10 / log(10 / Z0.LC); print('aeolian+neutral ustar is used')
    }  # reconstruct USTAR from U10 and Z0M
}
#plot.field(USTAR.aero.neu[,,8] / USTAR[,,8], lon.M2, lat.M2, type='def', zlim=c(0,2), col=TEMP_DIFF_65)
# then also reconstruct ustar under non-neutral condition
#plot.field(USTAR[,,6]*landfilt.M2, lon.M2, lat.M2, type='def', zlim=c(0,0.8), col=WBGYR)
#plot.field(USTAR.aeo.nonneu[,,6], lon.M2, lat.M2, type='def', zlim=c(0,0.8), col=WBGYR)
#plot.field(USTAR.aeo.nonneu[,,6]/USTAR.aeo.neu[,,6], lon.M2, lat.M2, type='def', zlim=c(0,2), col=TEMP_DIFF_65)
#plot.field(USTAR.aeo.neu[,,6], lon.M2, lat.M2, type='def', zlim=c(0,0.8), col=WBGYR)
#A = apply(USTAR, c(1,2), mean, na.rm=T)
#plot.field(A*landfilt.M2, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(0,0.7))
#A = apply(USTAR.aeo.nonneu, c(1,2), mean, na.rm=T)
#plot.field(A, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(0,0.7))
#A = apply(USTAR.aeo.neu, c(1,2), mean, na.rm=T)
#plot.field(A, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(0,0.7))

#print('finished threshold calculation; now calculate fluid and impact thresholds')

###################
# construct dust emission thresholds

# first get dry fluid threshold
ustar_ft.dry = 1.0*sqrt(0.0123 * (dns_slt*g*D_p + gamma/D_p)) / sqrt(RHOA)
#plot.field(ustar_ft.dry[,,6], lon.M2, lat.M2)

# to get wet fluid threshold,
# get gravimetric water content from volumetric water content
gwc_sfc = SHR_CONST_RHOFW*sweep(SFMC, c(1,2), bulk_den, '/')
#plot.field(gwc_sfc[,,6], lon.M2, lat.M2)
#plot.field.log(gwc_sfc[,,6], lon.M2, lat.M2, type='def',zlim=c(-3,0), legend.mar=5)

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
#plot.field(frc_thr_wet_fct[,,6], lon.M2, lat.M2, type='def', zlim=c(1,3.6))

# then get fluid and impact thresholds (now assume using FAO's moisture effect)
ustar_ft.wet = ustar_ft.dry * frc_thr_wet_fct
#ustar_it = ustar_ft.dry * B_it * landfilt.M2
ustar_it = B_it * sweep((ustar_ft.dry), c(1,2), landfilt.M2, '*')
#plot.field(ustar_it[,,6], lon.M2, lat.M2)
#plot.field(ustar_ft.wet[,,6], lon.M2, lat.M2)
#A = apply(USTAR - ustar_it, c(1,2), mean, na.rm=T)
#plot.field(A*landfilt.M2, lon.M2, lat.M2, legend.mar=5, type='sign', zlim=c(-0.5,0.5), col=TEMP_DIFF_65)
#A = apply(USTAR.aeo.nonneu - ustar_it, c(1,2), mean, na.rm=T)
#plot.field(A, lon.M2, lat.M2, legend.mar=5, type='sign', zlim=c(-0.5,0.5), col=TEMP_DIFF_65)

#print('finished threshold calculation; now calculate intermittency effect')

###################
# construct intermittency

# choose a friction velocity. Aeolian/Aerodynamic? Neutral/Non-neutral?
# USTAR.aeo.neu, USTAR.aeo.nonneu, USTAR.aero.neu, USTAR (aero+nonneu, default MERRA2 output)
for (aa in 1) {
  if (wind.scale=='default') {wnd_frc_slt = USTAR*F_eff.LC; print('default MERRA2 ustar is used for emission')}   #  default MERRA2 (aero+nonneu)output
  else if (wind.scale=='aeo+neu') {wnd_frc_slt = USTAR.aeo.neu*F_eff.LC; print('aeolian+neutral ustar is used for emission')}   #  reconstructed aeo+neu USTAR
  else if (wind.scale=='aeo+nonneu') {wnd_frc_slt = USTAR.aeo.nonneu*F_eff.LC; print('aeolian+neutral ustar is used for emission')}   #  reconstructed aeo+nonneu USTAR
  else if (wind.scale=='no.drag.partition') {wnd_frc_slt = USTAR; print('default MERRA2 ustar is used for emission but no drag partition effect')}
}

# if intermittency should be used, compute this section
if (use.intermittency=='TRUE') {
# get Monin-Obukhov length in m (Donald Golder, 1972; Khaled Essa, 2000)
# we use aerodynamic non-neutral USTAR here
Obu_L = - RHOA * c_p * T10M * USTAR^3 / (k * g * SHLAND)

u_mean_slt = (wnd_frc_slt/k) * log(0.1 / 1e-4)  # We used z0a = 1e-4 m in CESM and assumed. Should we change to aeolian data here?
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
#plot.field(u_mean_slt[,,6], lon.M2, lat.M2)
#plot.field(u_sd_slt[,,6], lon.M2, lat.M2)
#plot.field(u_fld_thr[,,6], lon.M2, lat.M2)
#plot.field(u_impct_thr[,,6], lon.M2, lat.M2)
#plot.field(intrmtncy_fct[,,6], lon.M2, lat.M2)
#A = apply(intrmtncy_fct, c(1,2), mean, na.rm=T)
#plot.field(A, lon.M2, lat.M2, col=WBGYR, legend.mar=5, type='def', zlim=c(0,1))

#print('finished intermittency; now calculate dust emission flux')

###################
# a few more essential variables for dust emission equations, modified 18 Nov 2021
if (DPM=='K14') {
# standardized fluid threshold (m/s)
ustar_st = ustar_ft.wet * sqrt(RHOA/rho_a0)
# dust emission coefficient (a measure of soil erodibility)
C_d = C_d0 * exp(-C_e * (ustar_st-ustar_st0)/ustar_st0);
#A = apply(C_d, c(1,2), mean, na.rm=T)
#plot.field.log(A/C_d0, lon.M2, lat.M2, type='def', zlim=c(-5,0), legend.mar=5)
# get fragmentation exponent and limit to frag_expt_lim and below
frag_expt = (C_alpha*(ustar_st-ustar_st0)/ustar_st0)
frag_expt[which(frag_expt>frag_expt_lim)] = frag_expt_lim
#A = apply(frag_expt, c(1,2), mean, na.rm=T)
#plot.field(A, lon.M2, lat.M2, legend.mar=5)
} #else if (DPM=='Z03') { # dml added 18 Nov 2021
# sandblasting efficiency
#dst_slt_flx_rat_ttl = 100*exp( log(10) * (13.4*mss_frc_cly_vld - 6))
#}

###################
# compute 3-D dust emissions

# compute dust emissions
for (aa in 1) {
if (DPM == 'K14') {
if (use.intermittency == 'TRUE') {
F_d.noeta = F_d = array(0, dim=dim(C_d))   # using dimension of C_d to define F_d dimension
ind = which(wnd_frc_slt > ustar_it)  # 3-D index for which wind speed > threshold
F_d.noeta[ind] = C_tune * C_d[ind] * f_bare[ind] * RHOA[ind] * (wnd_frc_slt[ind]^2-ustar_it[ind]^2)/ustar_it[ind] * (wnd_frc_slt[ind]/ustar_it[ind])^(frag_expt[ind])
F_d[ind] = F_d.noeta[ind] * intrmtncy_fct[ind]
#F_d[ind] = C_tune * C_d[ind] * f_bare[ind] * RHOA[ind] * (wnd_frc_slt[ind]^2-ustar_it[ind]^2)/ustar_it[ind] * (wnd_frc_slt[ind]/ustar_it[ind])^(frag_expt[ind]) * intrmtncy_fct[ind]             # calculate F_d for indices with wind speed > threshold
F_d = sweep(F_d, c(1,2), mss_frc_cly_vld, '*')  # multiply 3-D fields by 2-D fields
#F_d_sum = apply(F_d, c(1,2), sum, na.rm=T) * 3600 * dt * 365
#plot.field.log(F_d_sum, lon.M2, lat.M2, legend.mar=5, col=WBGYR)
} else if (use.intermittency == 'FALSE') {
F_d = array(0, dim=dim(C_d))   # using dimension of C_d to define F_d dimension
if (use.threshold == 'fluid') {
ind = which(wnd_frc_slt > ustar_ft.wet)  # 3-D index for which wind speed > threshold
print('use fluid threshold')
F_d[ind] = C_tune * C_d[ind] * f_bare[ind] * RHOA[ind] * (wnd_frc_slt[ind]^2-ustar_ft.wet[ind]^2)/ustar_ft.wet[ind] * (wnd_frc_slt[ind]/ustar_ft.wet[ind])^(frag_expt[ind])              # calculate F_d for indices with wind speed > threshold
} else if (use.threshold == 'impact') {
  ind = which(wnd_frc_slt > ustar_it)  # 3-D index for which wind speed > threshold
  F_d[ind] = C_tune * C_d[ind] * f_bare[ind] * RHOA[ind] * (wnd_frc_slt[ind]^2-ustar_it[ind]^2)/ustar_it[ind] * (wnd_frc_slt[ind]/ustar_it[ind])^(frag_expt[ind])              # calculate F_d for indices with wind speed > threshold
}

F_d = sweep(F_d, c(1,2), mss_frc_cly_vld, '*')  # multiply 3-D fields by 2-D fields
}

# save dust emission
#filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006_dt=3/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')
filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006_new/', directory, '/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')
filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006_new/', directory, '/intermittency_05x0625_', date_vec[d], '.RData', sep='')

for (aa in 1) {
if (use.intermittency == TRUE) {save('Obu_L', 'PBLH', 'intrmtncy_fct', 'ustar_ft.wet', 'ustar_it', 'wnd_frc_slt', 'lon.M2', 'lat.M2', file=filename)}#; save('F_d', 'F_d.noeta', 'intrmtncy_fct', 'lon.M2', 'lat.M2', file=filename)} 
  else save('F_d', 'lon.M2', 'lat.M2', file=filename)
}

# save thresholds
#filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/thresholds_05x0625_', date_vec[d], '.RData', sep='')
#save('ustar_ft.wet', 'ustar_it', 'lon.M2', 'lat.M2', file=filename)

}
else if (DPM == 'Z03') { # dml added 18 Nov 2021
wnd_frc_rat = ustar_ft.wet / wnd_frc_slt
F_d = array(0, dim=dim(ustar_ft.wet))   # using dimension of ustar_ft.wet to define F_d dimension
ind = which(wnd_frc_slt > ustar_ft.wet)  # 3-D index for which wind speed > threshold
F_d[ind] = flx_mss_fdg_fct * cst_slt * f_bare[ind] * RHOA[ind] * (wnd_frc_slt[ind]^3) * (1 - wnd_frc_rat[ind]) * (1 + wnd_frc_rat[ind]) * (1 + wnd_frc_rat[ind]) / g
F_d = sweep(F_d, c(1,2), dst_slt_flx_rat_ttl*mbl_bsn_fct , '*')  # multiply 3-D fields by 2-D fields
filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_',DPM,'_src_fn=',src_fn,'/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')
save('F_d', 'lon.M2', 'lat.M2', file=filename)
print('computed Z03 emission')
}

}

print('emission data saved')

}
proc.time() - ptm
#A = apply(F_d, c(1,2), sum, na.rm=T) *dt*3600*365
#plot.field.log(A ,lon.M2, lat.M2)



###################################################
###################################################
# one small loop for computing and saving hybrid drag partition factor F_eff
for (d in length(date_vec)) {# d-th day
  print(paste('day ', d, ': ', date_vec[d], sep=''))
  # get 3D data in time and space, subsetted as regional data with variable time resolutions
  i1 = indi[1]; j1 = indj[1]; t1 = indt[1]
  il = length(indi); jl = length(indj); tl = length(indt)
  #ptm = proc.time()
  nc = nc_open(files_lnd[d])#; names(nc$var)
  LAI = ncvar_get(nc, 'LAI', start=c(i1,j1,1), count=c(il,jl,24))[,,6]#[,,indt]
  # get drag partition factors and hybrid roughness length (later on, make this as a separate time loop independent of dust emission calculation?)
  f_bare = 1 - LAI/LAI_thr; f_bare[which(f_bare<0)] = NaN
  # Okin's vegetation drag partition
  K = pi/2 * (1/LAI - 1); K[which(K<0)] = 0
  SSR = (K + U_0 * c) / (K + c)  # Pierre-Okin shear stress ratio (SSR)
  # combine rock and vegetation drag partition factors
  if (LC=='LC1') {
    f_eff_veg.comp = sweep(SSR^3, c(1,2), (frc.veg.M2), '*')
    f_eff_mix.comp = sweep(SSR^3, c(1,2), (frc.mix.M2*(f_eff.r^3)), '*')
    F_eff.LC = (sweep(f_eff_veg.comp+f_eff_mix.comp, c(1,2), (frc.r.M2*f_eff.r^3), '+')) ^ (1/3)
    rm('f_eff_veg.comp','f_eff_mix.comp')  
  }
  else if (LC=='LC0') {
  f_eff_veg.comp = sweep(SSR^3, c(1,2), (frc.veg.M2+frc.mix.M2), '*')
  F_eff.LC = (sweep(f_eff_veg.comp, c(1,2), (frc.r.M2*f_eff.r^3), '+')) ^ (1/3)
  rm(f_eff_veg.comp)
  }
  #F_eff.LC = (frc.r.M2*f_eff.r^3 +  )^(1/3) 
  # save daily data (diurnal fluctuations are too small)
  #filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_F_eff_',LC,'/F_eff_hybrid_05x0625_', date_vec[d], '.RData', sep='')
  filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006_new/', directory, '/F_eff_hybrid_05x0625_2006_', date_vec[d], '.RData', sep='')
  save('F_eff.LC', 'SSR', 'LAI', 'lon.M2', 'lat.M2', file=filename)
}

# plot the required variables
.env = new.env()
F_eff.LC = SSR = LAI = array(0, dim=c(length(lon.M2), length(lat.M2)))
for (d in 1:length(date_vec)) {# d-th day
  print(date_vec[d])
  #filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_F_eff_',LC,'/F_eff_hybrid_05x0625_', date_vec[d], '.RData', sep='')
  filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006_new/', directory, '/F_eff_hybrid_05x0625_2006_', date_vec[d], '.RData', sep='')
  load(filename, envir = .env)
  F_eff.LC = F_eff.LC + .env$F_eff.LC
  SSR = SSR + .env$SSR
  LAI = LAI + .env$LAI
}
SSR = SSR/365; LAI = LAI/365; F_eff.LC = F_eff.LC/365

i.M2 = c(20:576,1:19)
lon.M2.r = lon.M2[i.M2]; lon.M2.r[558:576]=lon.M2.r[558:576]+360
LAI.1 = LAI; LAI.1[which(LAI>1)]=NaN
plot.field(LAI.1[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(0,1))
SSR.1 = SSR; SSR.1[which(LAI>1)]=NaN
plot.field(SSR.1[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(0,1))
plot.field(SSR.1[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(0.3,1))

F_eff.LC.1 = F_eff.LC; F_eff.LC.1[which(LAI>1)] = NaN
plot.field(F_eff.LC.1[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(0.3,1), col=WBGYR)
plot.field(F_eff.LC.1[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(0.3,1))

plot.field.proj(LAI[i.M2,], lon.M2.r, lat.M2, def.zlim=T, zlim=c(0,1))
plot.field.proj(SSR[i.M2,], lon.M2.r, lat.M2, def.zlim=T, zlim=c(0.3,1))
plot.field.proj(F_eff.LC[i.M2,], lon.M2.r, lat.M2, def.zlim=T, zlim=c(0,1))
plot.field.proj(F_eff.LC[i.M2,], lon.M2.r, lat.M2, def.zlim=T, zlim=c(0,1), col=WBGYR)

filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006_new/', directory, '/F_eff_hybrid_05x0625_2006_annual_mean.RData', sep='')
save('SSR', 'SSR.1', 'LAI', 'LAI.1', 'F_eff.LC', 'F_eff.LC.1' , 'lon.M2', 'lat.M2', file=filename)

# save daily dust emission thresholds
ptm = proc.time()
#for (d in 1) {# d-th day
for (d in 1:length(date_vec)) {# d-th day
  
  print(paste('day ', d, ': ', date_vec[d], sep=''))
  #print('extract MERRA2 data', sep='')
  ###################
  # get 3D data in time and space, subsetted as regional data with variable time resolutions
  i1 = indi[1]; j1 = indj[1]; #t1 = indt[1]
  il = length(indi); jl = length(indj); #tl = length(indt)
  
  
  #ptm = proc.time()
  nc = nc_open(files_flx[d])#; names(nc$var)#, start=c(11,26), count=c(5,5)
  #USTAR = ncvar_get(nc, 'USTAR', start=c(i1,j1,1), count=c(il,jl,24)); USTAR = USTAR[,,indt] # m/s
  #Z0M = ncvar_get(nc, 'Z0M', start=c(i1,j1,1), count=c(il,jl,24)); Z0M = Z0M[,,indt]     # m
  RHOA = ncvar_get(nc, 'RHOA', start=c(i1,j1,1), count=c(il,jl,24)); RHOA = RHOA[,,indt]     # kg/m3
  #PBLH = ncvar_get(nc, 'PBLH', start=c(i1,j1,1), count=c(il,jl,24)); PBLH = PBLH[,,indt]     # kg/m3
  nc = nc_open(files_lnd[d])#; names(nc$var)
  SFMC = ncvar_get(nc, 'SFMC', start=c(i1,j1,1), count=c(il,jl,24)); SFMC = SFMC[,,indt]
  
  ###################
  # construct dust emission thresholds
  
  # first get dry fluid threshold
  ustar_ft.dry = 1.0*sqrt(0.0123 * (dns_slt*g*D_p + gamma/D_p)) / sqrt(RHOA)
  #plot.field(ustar_ft.dry[,,6], lon.M2, lat.M2)
  
  # to get wet fluid threshold,
  # get gravimetric water content from volumetric water content
  gwc_sfc = SHR_CONST_RHOFW*sweep(SFMC, c(1,2), bulk_den, '/')
  #plot.field(gwc_sfc[,,6], lon.M2, lat.M2)
  #plot.field.log(gwc_sfc[,,6], lon.M2, lat.M2, type='def',zlim=c(-3,0), legend.mar=5)
  
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
  #plot.field(frc_thr_wet_fct[,,6], lon.M2, lat.M2, type='def', zlim=c(1,3.6))
  
  # then get fluid and impact thresholds (now assume using FAO's moisture effect)
  ustar_ft.wet = ustar_ft.dry * frc_thr_wet_fct
  #ustar_it = ustar_ft.dry * B_it * landfilt.M2
  ustar_it = B_it * sweep((ustar_ft.dry), c(1,2), landfilt.M2, '*')
 
  # save thresholds
  filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006_new/', directory, '/thresholds_05x0625_', date_vec[d], '.RData', sep='')
  save('ustar_it', 'ustar_ft.wet', 'lon.M2', 'lat.M2', file=filename) 
}