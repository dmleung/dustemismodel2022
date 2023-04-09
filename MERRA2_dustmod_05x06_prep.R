
# MERRA2 processing trial
# 25 Sep 2021
# Main script for the first project

setwd('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR')
# libraries (optional?)
library(ncdf4); library(fields); library(maps); library(Metrics); library(pracma)
# functions
source("get_geo.R"); source("get_met.R"); source("sptial_plot_fns.R")
source("dust_research_fns.R")
# color scheme
load("WBGYR_scheme.RData"); load('TEMP_DIFF.RData')


#################################################################
n = 18  # specify timestep

#nc = nc_open('/Volumes/GoogleDrive/My Drive/MERRA2/upscale_2006/tavg1_2d_flx_Nx/MERRA2_300.tavg1_2d_flx_Nx.20060101.SUB.nc')
nc = nc_open('/Volumes/GoogleDrive/My Drive//MERRA2/upscale_2006/MERRA2_300.tavg1_2d_lnd_Nx.20060101.SUB.nc')
lon.M2 = ncvar_get(nc, 'lon'); lat.M2 = ncvar_get(nc, 'lat')
nc_close(nc)

# 2D lat and lon
lon_mat = matrix(lon.M2, nrow = length(lon.M2),ncol = length(lat.M2)); lat_mat = t(matrix(lat.M2, nrow = length(lat.M2),ncol = length(lon.M2)))


names(nc$var)
USTAR = ncvar_get(nc, 'USTAR')[,,n]
Z0M = ncvar_get(nc, 'Z0M')[,,n]
RHOA = ncvar_get(nc, 'RHOA')[,,n]
dim(USTAR)
nc = nc_open('/Volumes/SEAGATE/MERRA2/upscale_2006/tavg1_2d_slv_Nx/MERRA2_300.tavg1_2d_slv_Nx.20060101.SUB.nc')
names(nc$var)
DISPH = ncvar_get(nc, 'DISPH')[,,n]
T10M = ncvar_get(nc, 'T10M')[,,n]
U10M = ncvar_get(nc, 'U10M')[,,n]
V10M = ncvar_get(nc, 'V10M')[,,n]
U10 = sqrt(U10M^2 + V10M^2)  # define U10 as the total wind speed of U10M and V10M. Don't confuse U10 with U10M.

plot.field(T10M, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(220,320))

#define spatial domain
indi = 251:390; indj = 171:260   
plot.field(RHOA[indi,indj], lon.M2[indi], lat.M2[indj])
plot.field(DISPH[indi,indj], lon.M2[indi], lat.M2[indj])
plot.field(USTAR[indi,indj], lon.M2[indi], lat.M2[indj])
plot.field(U10M[indi,indj], lon.M2[indi], lat.M2[indj])
plot.field(V10M[indi,indj], lon.M2[indi], lat.M2[indj])
plot.field(U10[indi,indj], lon.M2[indi], lat.M2[indj])
plot.field(Z0M[indi,indj], lon.M2[indi], lat.M2[indj])
plot.map.log(Z0M[indi,indj],lon.M2[indi], lat.M2[indj],type='def',zlim=c(-3,-0), legend.mar=5)

plot.map.log((Z0M[indi,indj])^2,lon.M2[indi], lat.M2[indj],type='def',zlim=c(-10,0), legend.mar=5)
plot.map.log((USTAR*Z0M)^2,lon.M2, lat.M2,type='def',zlim=c(-10,1), legend.mar=5)
plot.map.log(DISPH,lon.M2, lat.M2,type='def',zlim=c(-3,2), legend.mar=5)
plot.field(DISPH, lon.M2, lat.M2)
plot.field(DISPH[indi,indj], lon.M2[indi], lat.M2[indj])
# use displacement height as landfilt (same for other reoslution; use displacement height as landfilt)
landfilt.M2 = DISPH; landfilt.M2 = landfilt.M2/landfilt.M2
image(landfilt.M2)
quantile(landfilt.M2,na.rm=T)
plot.map.log(Z0M*landfilt.M2*100, lon.M2, lat.M2, type='def',zlim=c(-3,3), legend.mar=5)
save('landfilt.M2','lon.M2','lat.M2',file='landfilt_05x0625.RData')

# test if U10, USTAR, and Z0M are consistent with each other
USTAR.aero.neu = 0.4 * U10 / log(10 / Z0M)  # reconstruct USTAR from U10 and Z0M
plot.field(USTAR.aero.neu[indi,indj], lon.M2[indi], lat.M2[indj], type='def', zlim=c(0,1))
plot.field(USTAR[indi,indj], lon.M2[indi], lat.M2[indj], type='def', zlim=c(0,1))
plot.field(USTAR.aero.neu[indi,indj] / USTAR[indi,indj], lon.M2[indi], lat.M2[indj], type='def', zlim=c(0,2), col=TEMP_DIFF_65)
plot.map.log(USTAR.aero.neu[indi,indj] / USTAR[indi,indj], lon.M2[indi], lat.M2[indj], type='def',zlim=c(-1,1), legend.mar=5, col=TEMP_DIFF_65)


# comparing aerodynamic ustar with aeolian ustar
.env = new.env()
load('/Volumes/GoogleDrive/My Drive/Aeolianz0ErodDust4astrid/Pr05_z0.tran_199701-199712_025x025.RData', envir=.env)
lon.vec = .env$lon.vec; lat.vec = .env$lat.vec
z0.tran = .env$z0.tran; 
indi.025 = 470:820+160; indj.025 = 340:520
z0.tran = z0.tran[indi.025,indj.025,]
lon.vec = lon.vec[indi.025]; lat.vec = lat.vec[indj.025]
dim(z0.tran)
plot.map.log(z0.tran[,,1]*100, lon.vec, lat.vec, type='def', zlim=c(-3,3), legend.mar=5)

# regrid aeolian roughness first
z0.05_1 = sp.dissolve.fast1(z0.tran[,,1], lon.vec, lat.vec, lon.M2[indi], lat.M2[indj])
plot.map.log(z0.05_1*landfilt.M2[indi,indj]*100, lon.M2[indi], lat.M2[indj], type='def', zlim=c(-3,3), legend.mar=5)

# compare Pr05 aeolian roughness with MERRA2 aerodynamic roughness
plot.map.log(Z0M[indi,indj]*100, lon.M2[indi], lat.M2[indj],type='def',zlim=c(-3,1), legend.mar=5)
plot.map.log(z0.05_1, lon.M2[indi], lat.M2[indj], type='def',zlim=c(-3,1), legend.mar=5)


# calculate ustar
KSI_m = 0.4*U10/USTAR - log(10 / Z0M) # first determine MOST instability (KSI)
USTAR.aeo.nonneu = 0.4 * U10[indi,indj] / (log(1000 / z0.05_1)+KSI_m[indi,indj])
USTAR.aeo.neu = 0.4 * U10[indi,indj] / log(1000 / z0.05_1)  # reconstruct USTAR from U10 and Z0M
plot.field(USTAR.aero.neu[indi,indj]*landfilt[indi,indj], lon.M2[indi], lat.M2[indj], type='def', zlim=c(0,1))
plot.field(USTAR.aeo.neu, lon.M2[indi], lat.M2[indj], type='def', zlim=c(0,1))
plot.field(USTAR.aeo.nonneu, lon.M2[indi], lat.M2[indj], type='def', zlim=c(0,1))

plot.map.log(z0.05_1 / (Z0M[indi,indj]*100), lon.M2[indi], lat.M2[indj], type='def',zlim=c(-3,0), legend.mar=5)
#plot.map.log(USTAR.aeo / (USTAR.rec[indi,indj]), lon.M2[indi], lat.M2[indj], type='def',zlim=c(-0.307,-0.045), legend.mar=5, ticks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
plot.field(USTAR.aeo / (USTAR.rec[indi,indj]), lon.M2[indi], lat.M2[indj], type='def', zlim=c(0.5,0.9))

######################################################################################
# constant data
# land fraction, ocean, landice, and lake fractions, and surface geopotential height

nc = nc_open('/Volumes/GoogleDrive/My Drive/MERRA2/upscale_2006/MERRA2_101.const_2d_asm_Nx.00000000.nc4')
names(nc$var)
FROCEAN = ncvar_get(nc, 'FROCEAN') # fraction ocean
FRLANDICE = ncvar_get(nc, 'FRLANDICE') # fraction land ice
PHIS = ncvar_get(nc, 'PHIS') # geopotential
FRLAKE = ncvar_get(nc, 'FRLAKE') # fraction lake
FRLAND = ncvar_get(nc, 'FRLAND') # fraction land 

plot.field(PHIS[indi,indj]/9.81, lon.M2[indi], lat.M2[indj], legend.mar = 4)
plot.field(PHIS[iC05-288,90:170], lon.M2[iC05-288], lat.M2[90:170], legend.mar=5)
#Eplot.field(PHIS/9.81, lon.M2, lat.M2, legend.mar = 4)
plot.field(FRLANDICE, lon.M2, lat.M2)
plot.field(FRLAND[indi,indj], lon.M2[indi], lat.M2[indj])


################################################################
# calculate dry fluid threshold using RHOA (in m s^-1)
g = 9.81        # m/s^2, gravitational acceleration 
dns_slt = 2650  # kg/m3, soil particle density
D_p = 130e-6    # m, soil median diameter
ustar_ft.dry = 1.0*sqrt(0.0123 * (dns_slt*g*D_p + 5.0e-4/D_p)) / sqrt(RHOA)
plot.field(ustar_ft.dry[indi,indj], lon.M2[indi], lat.M2[indj])

# get porosity
nc = nc_open('/Volumes/GoogleDrive/My Drive/MERRA2/upscale_2006/MERRA2_100.const_2d_lnd_Nx.00000000.nc4')
names(nc$var)
poros = ncvar_get(nc, 'poros') # porosity (m3 / m3, also known as volumetric soil water at saturation)
plot.field(poros[indi,indj], lon.M2[indi], lat.M2[indj])
# get bulk density of dry soil material [kg/m^3]
bulk_den = (1 - poros)*dns_slt 
plot.field(bulk_den[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5)
# get land-surface variables
nc = nc_open('/Volumes/SEAGATE/MERRA2/upscale_2006/tavg1_2d_lnd_Nx/MERRA2_300.tavg1_2d_lnd_Nx.20060101.SUB.nc')
names(nc$var)
SFMC = ncvar_get(nc, 'SFMC')[,,n] # volumetric water content (m3 / m3)
LAI = ncvar_get(nc, 'LAI')[,,n] # LAI (m2 / m2)
SHLAND = ncvar_get(nc, 'SHLAND')[,,n] # LAI (W / m2)
plot.field(LAI[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5)
plot.field(LAI, lon.M2, lat.M2, legend.mar=5)
plot.field(SHLAND[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5)
plot.field(SHLAND, lon.M2, lat.M2, legend.mar=5, type='sign', col=TEMP_DIFF_65, zlim=c(-500,500))
plot.field(SFMC[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5)
plot.field(SFMC, lon.M2, lat.M2)
plot.field.log(SFMC, lon.M2, lat.M2, legend.mar=5)
# get gravimetric water content from volumetric water content
SHR_CONST_RHOFW = 1000   # kg / m3, pure water density
gwc_sfc = SFMC * SHR_CONST_RHOFW / bulk_den   # surface gravimetric water content (kg/kg)
plot.field(gwc_sfc[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5)
plot.field(gwc_sfc, lon.M2, lat.M2)
plot.map.log(gwc_sfc, lon.M2, lat.M2, type='def',zlim=c(-3,0), legend.mar=5)

# check if MERRA2 and CESM soil moisture are consistent with each other
nc = nc_open('CESMmetinput/3hourly/extractmetforR.clm2.h3.2012-01-01-00000.nc')
lon_19x25 = seq(-180,177.5, by=2.5)
lat_19x25 = ncvar_get(nc, 'lat')
names(nc$var)
GWC = ncvar_get(nc, 'GWC')
ntime = dim(GWC)[3]  # ntime = 8, for 3-hourly data
GWC = GWC[c(73:144,1:72),,1]
plot.field(GWC[60:100,40:75], LON[60:100], LAT[40:75])
plot.field(GWC, LON, LAT)
plot.map.log(GWC, LON, LAT, type='def',zlim=c(-3,0), legend.mar=5)
# use GWC to get a landfilt for 1.9x2.5 CESM
landfilt.CESM.19x25 = GWC/GWC
#image(landfilt.CESM.19x25)

# then get moisture factor
# use SoilGrids' or CESM's FAO clay fraction
nc = nc_open('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR/CESMmetinput/surfdata_1.9x2.5_simyr1850_c091108.nc')  # CESM's FAO (2012) clay
f_clay.FAO.19x25 = ncvar_get(nc, 'PCT_CLAY')[c(73:144,1:72),,1]*0.01  # 0-1, fraction of clay at 1st (surface) soil level
f_clay.FAO.19x25 = f_clay.FAO.19x25*landfilt.CESM.19x25 # filter out ocean filled values
plot.field(f_clay.FAO.19x25, LON, LAT)
#title('clay fraction at 1st nlevsoi: surface')
# bilinear downscale f_clay.FAO.19x25 to MERRA2 resolution
f_clay.FAO.05x06 = sp.smooth(f_clay.FAO.19x25, LON, LAT, lon.M2, lat.M2, landfilt=landfilt.M2)
plot.field(f_clay.FAO.05x06, lon.M2, lat.M2, type='def', zlim=c(0,0.65))
save('f_clay.FAO.05x06', file='f_clay_FAO_05x06.RData')

# compare with SoilGrids soil texture
#head(load('/Volumes/GoogleDrive/My Drive/soil_database/SoilGrids_texture.RData'))
load('/Volumes/GoogleDrive/My Drive/soil_database/SoilGrids_texture.RData', envir=.env)
ls(.env)
lon.SG = .env$lon; lat.SG = .env$lat   # about 0.2502 deg lon and 0.2505 deg lat
f_clay.SG.025 = (.env$texture)[,,7]*0.01
f_clay.SG.025 = fill.miss.data(f_clay.SG.025,surrounding=5)
plot.field(f_clay.SG.025, lon.SG, lat.SG, type='def', zlim=c(0,0.65))
f_clay.SG.05x06 = sp.dissolve.fast(f_clay.SG.025, lon.SG, lat.SG, lon.M2, lat.M2)
plot.field(f_clay.SG.05x06, lon.M2, lat.M2, type='def', zlim=c(0,0.65))
save('f_clay.SG.05x06', file='f_clay_SoilGrids_05x06.RData')
# threshold gravimetric moisture and clay limitation on dust emission
gwc_thr.SG = 0.01*(0.17*f_clay.SG.05x06*100 + 0.0014*(f_clay.SG.05x06*100)^2)#Fecan et al. (1999)
mss_frc_cly_vld.SG = f_clay.SG.05x06; mss_frc_cly_vld.SG[which(mss_frc_cly_vld.SG>0.2)] = 0.2
plot.field(gwc_thr.SG, lon.M2, lat.M2)

load('f_clay_FAO_05x06.RData')
gwc_thr.FAO = 0.01*(0.17*f_clay.FAO.05x06*100 + 0.0014*(f_clay.FAO.05x06*100)^2)
mss_frc_cly_vld.FAO = f_clay.FAO.05x06; mss_frc_cly_vld.FAO[which(mss_frc_cly_vld.FAO>0.2)] = 0.2
plot.field(gwc_thr.FAO, lon.M2, lat.M2)
plot.field(mss_frc_cly_vld.FAO, lon.M2, lat.M2)

# moisture effect on fluid threshold
frc_thr_wet_fct.SG = sqrt(1.0 + 1.21*(100.0*(gwc_sfc - gwc_thr.SG))^0.68); frc_thr_wet_fct.SG[which(gwc_sfc<gwc_thr.SG)] = 1

frc_thr_wet_fct.FAO = sqrt(1.0 + 1.21*(100.0*(gwc_sfc - gwc_thr.FAO))^0.68); frc_thr_wet_fct.FAO[which(gwc_sfc<gwc_thr.FAO)] = 1
plot.field(frc_thr_wet_fct.FAO, lon.M2, lat.M2)
plot.field(frc_thr_wet_fct.SG, lon.M2, lat.M2)

# then get fluid and impact thresholds (now assume using FAO's moisture effect)
B_it = 0.81
ustar_ft.wet = ustar_ft.dry * frc_thr_wet_fct.FAO
ustar_it = ustar_ft.dry * B_it * landfilt.M2
plot.field(ustar_it, lon.M2, lat.M2)
plot.field(ustar_ft.wet, lon.M2, lat.M2)

#############################################################
# use LAI to get bare land fraction using Natalie Mahowald's equation
# and vegetation drag partition using 
nc = nc_open('/Volumes/SEAGATE/MERRA2/upscale_2006/tavg1_2d_lnd_Nx/MERRA2_300.tavg1_2d_lnd_Nx.20060101.SUB.nc')
names(nc$var)
LAI = ncvar_get(nc, 'LAI')[,,n] # LAI (m2 / m2)
plot.field(LAI[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5)
plot.field(LAI, lon.M2, lat.M2)
plot.map.log(LAI, lon.M2, lat.M2, legend.mar=5, zlim=c(-2,1))

# Natalie's bare land fraction approximation
LAI_thr = 1   # Mahowald et al. (2010) set as 0.3 but we changed to 1
f_bare = 1 - LAI/LAI_thr; f_bare[which(f_bare<0)] = NaN
plot.field(f_bare, lon.M2, lat.M2)
# make a filter for dust emission sources
emisfilt = f_bare/f_bare    # note that this filt changes with time, or one should make a annual one
#image(emisfilt)


# Okin's vegetation drag partition
U_0 = 0.32  # SSR in the immediate lee of a plant, dimensionless
c = 4.8     # e-folding distance SSR recovery, dimensionless
# Using Caroline Pierre (2014)'s formulation of Okin's scheme
K = pi/2 * (1/LAI - 1); K[which(K<0)] = 0
SSR = (K + U_0 * c) / (K + c)  # Pierre-Okin shear stress ratio (SSR)
plot.field(SSR, lon.M2, lat.M2)

# get rock drag partition from Pr05 min z0 data
load('/Volumes/GoogleDrive/My Drive/Aeolianz0ErodDust4astrid/Pr05_z0.tran_199701-199712_025x025.RData', envir=.env)
ls(.env)
lon.025 = .env$lon.vec; lat.025 = .env$lat.vec
z0.tran = .env$z0.tran   # data in cm
z0.min = apply(z0.tran, c(1,2), min, na.rm=T); z0.min[which(is.infinite(z0.min))]=NaN
z0.M2 = sp.dissolve.fast(z0.min, lon.025, lat.025, lon.M2, lat.M2)
#plot.map.log(z0.min, lon.025, lat.025, legend.mar=5)
plot.map.log(z0.M2, lon.M2, lat.M2, type='def', zlim=c(-3,0), col=tim.colors(32)[32:1], legend.mar=5)
rm(z0.tran)   # save memory
save('z0.M2', file='Pr05_z0.min_05x0625.RData')

load('Pr05_z0.min_05x0625.RData')
# in this formulation, set globally constant z0s because soil median diameter ~ 130 um
z0s = D_p / 15    # m, smooth roughness length, White (2006) equation
X = 10            # m, we assume interobstacle distance ~ 10 m
f_eff.r = 1 - ( log(z0.M2/100/z0s) / log(0.7*(10/z0s)^0.8) )  
# limit f_eff in between 0 and 1
#f_eff.r[which(f_eff.r>0.65 & lon_mat>80)] = f_eff.r[which(f_eff.r>0.65 & lon_mat>80)] + 0.1
f_eff.r[which(f_eff.r>=1)] = 1
f_eff.r[which(f_eff.r<=0)] = 0
plot.field(f_eff.r, lon.M2, lat.M2, type='def', zlim=c(0.3,1))
#i = 200:400; j = 160:255 # N Africa + Middle East

i = 236:345; j = 110:255 
#i = 236:400; j = 130:255
#i = 236:400; j = 160:255
f_eff.r[i,j] = f_eff.r[i,j]*0.95
i = 346:400; j = 110:210
f_eff.r[i,j][which(f_eff.r[i,j]>0.75)]=0.01
f_eff.r[i,j] = f_eff.r[i,j]*0.95
#i = 410:510; j = 255:290
#i = 410:510; j = 255:310
i = 410:510; j = 255:275
#i = 370:520; j = 215:320
#plot.field(f_eff.r[i,j], lon.M2[i], lat.M2[j], type='def', zlim=c(0.3,1))
f_eff.r[i,j] = f_eff.r[i,j]*1.05
i = 430:510; j = 255:310
f_eff.r[i,j] = f_eff.r[i,j]*1.05
#i = 425:510; j = 255:275
#f_eff.r[i,j][which(f_eff.r[i,j]<0.75)]=f_eff.r[i,j][which(f_eff.r[i,j]<0.75)]*1.05

plot.field(f_eff.r, lon.M2, lat.M2, type='def', zlim=c(0.3,1))
plot.field(f_eff.r[i,j], lon.M2[i], lat.M2[j], type='def', zlim=c(0.3,1))
plot.field.log(z0.M2, lon.M2, lat.M2, type='def', zlim=c(-3,1))
# compare with Pr12 data
head(load('f_eff_05x05_globe_SoilGrids_Dp.RData'))
load('f_eff_05x05_globe_SoilGrids_Dp.RData',envir=.env)
lon.05 = .env$lon_re; lat.05 = .env$lat_re
f_eff.new.M2 = sp.dissolve(f_eff.new, lon.05, lat.05, lon.M2, lat.M2)
plot.field(f_eff.new.M2, lon.M2, lat.M2, type='def', zlim=c(0,1))
plot.field(f_eff/f_eff.new.M2, lon.M2, lat.M2, type='def', zlim=c(0.5,1.5), col=TEMP_DIFF_65)  # it seems the two drag partition maps are very close to each other.


# Then load GLCNMO data and combine with drag partition factors using LC equations.
# when regridding, only need to do it once. So do it first globally, and then save globally and later extract globally or regionally when needed.
head(load('GLCNMO_frc_area.RData'))
load('GLCNMO_frc_area.RData',envir=.env)
ls(.env)
lon.rf = .env$lon.rf; lat.rf = .env$lat.rf; frc = .env$frc
frc.veg = apply(frc[,,c(6:8,11:13)], c(1,2), sum, na.rm=T)
frc.veg.M2 = sp.dissolve(frc.veg, lon.rf, lat.rf, lon.M2, lat.M2)
#frc.tree = apply(frc[,,6:7], c(1,2), sum, na.rm=T)
frc.r = apply(frc[,,c(16:17)], c(1,2), sum, na.rm=T); 
frc.r.M2 = sp.dissolve(frc.r, lon.rf, lat.rf, lon.M2, lat.M2)
frc.mix = apply(frc[,,c(9:10)], c(1,2), sum, na.rm=T)
frc.mix.M2 = sp.dissolve(frc.mix, lon.rf, lat.rf, lon.M2, lat.M2)
rm(frc.mix, frc.r, frc.veg, frc, lon.rf, lat.rf)   # save space and memory
plot.field(frc.mix.M2, lon.M2, lat.M2, type='def', zlim=c(0,1), col=WBGYR)
plot.field(frc.r.M2, lon.M2, lat.M2, type='def', zlim=c(0,1), col=WBGYR)
plot.field(frc.veg.M2, lon.M2, lat.M2, type='def', zlim=c(0,1), col=WBGYR)
plot.field(frc.veg.M2+frc.mix.M2, lon.M2, lat.M2, type='def', zlim=c(0,1), col=WBGYR)
plot.field((frc.veg.M2+frc.r.M2+frc.mix.M2)*emisfilt, lon.M2, lat.M2, type='def', zlim=c(0,1), col=WBGYR)
save('frc.r.M2', 'frc.mix.M2', 'frc.veg.M2', 'lon.M2', 'lat.M2', file='GLCNMO_frc_area_05x0625.RData')

load('GLCNMO_frc_area_05x0625.RData')
# then combine drag partition factors
F_eff.LC1 = (frc.r.M2*f_eff.r^3 + frc.veg.M2*SSR^3 + frc.mix.M2*(SSR*f_eff.r)^3)^(1/3)
F_eff.LC0 = (frc.r.M2*f_eff.r^3 + (frc.veg.M2+frc.mix.M2)*SSR^3)^(1/3)  # newly defined LC0, removed mixed regime
plot.field(F_eff.LC1, lon.M2, lat.M2, type='def', zlim=c(0,1), col=WBGYR)
plot.field(F_eff.LC0, lon.M2, lat.M2, type='def', zlim=c(0,1))

# construct hybrid roughness length (m) and hybrid aeolian ustar
Z0.LC2 = z0s * exp( (1 - F_eff.LC0) * log(0.7*(X/z0s)^0.8) )   # in m
plot.field.log(Z0.LC2,lon.M2,lat.M2,legend.mar=5,col=tim.colors(32)[32:1])
plot.field.log(Z0M*landfilt.M2,lon.M2,lat.M2,legend.mar=5,col=tim.colors(32)[32:1])
KSI_m = 0.4*U10/USTAR - log(10 / Z0M) # first determine MOST instability (KSI)
USTAR.aeo.nonneu = 0.4 * U10 / (log(1000 / Z0.LC2)+KSI_m)
USTAR.aeo.neu = 0.4 * U10 / log(1000 / Z0.LC2)  # reconstruct USTAR from U10 and Z0.LC
#plot.field(USTAR.aero.neu[indi,indj]*landfilt[indi,indj], lon.M2[indi], lat.M2[indj], type='def', zlim=c(0,1))
plot.field(USTAR.aeo.neu[indi,indj], lon.M2[indi], lat.M2[indj], type='def', zlim=c(0,1))
plot.field(USTAR.aeo.nonneu, lon.M2, lat.M2, type='def', zlim=c(0,1))

plot.map.log(Z0.LC2 / (Z0M*100), lon.M2, lat.M2, type='def',zlim=c(-5,-2), legend.mar=5)
#plot.map.log(USTAR.aeo / (USTAR.rec[indi,indj]), lon.M2[indi], lat.M2[indj], type='def',zlim=c(-0.307,-0.045), legend.mar=5, ticks=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
plot.field(USTAR.aeo.nonneu / (USTAR), lon.M2, lat.M2, type='def', zlim=c(0,1))

plot.field(USTAR.aeo.nonneu, lon.M2, lat.M2, type='def', zlim=c(0,1))
plot.field(USTAR*landfilt.M2, lon.M2, lat.M2, type='def', zlim=c(0,1))

plot.field(USTAR.aeo.nonneu / (USTAR) * F_eff.LC0, lon.M2, lat.M2, type='def', zlim=c(0,1))
plot.field(KSI_m*landfilt.M2, lon.M2, lat.M2, type='def', zlim=c(0,20))
plot.field(KSI_m*landfilt.M2, lon.M2, lat.M2, type='sign', zlim=c(-20,20))
plot.field(log(1000 / Z0.LC2), lon.M2, lat.M2, type='def', zlim=c(0,20))

#############################################################
# For Z03, need to process source functions (soil erodibility maps) too
# added 18 Nov 2021
# Ginoux's source function
nc = nc_open('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR/ginoux_dust_source_0.25x0.25.nc'); names(nc$var)
lon.G025 = ncvar_get(nc, 'lon'); lat.G025 = ncvar_get(nc, 'lat')
mbl_bsn_fct_G.025x025 = ncvar_get(nc,'source')
plot.field(mbl_bsn_fct_G.025x025, lon.G025, lat.G025)
mbl_bsn_fct_G.05x06 = sp.dissolve(mbl_bsn_fct_G.025x025, lon.G025, lat.G025, lon.M2, lat.M2)
plot.field(mbl_bsn_fct_G.05x06, lon.M2, lat.M2)
plot.field.log(mbl_bsn_fct_G.05x06, lon.M2, lat.M2)
save('mbl_bsn_fct_G.05x06', 'lon.M2', 'lat.M2', file='ginoux_dust_source_05x0625.RData')
# Zender's source function
nc = nc_open('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR/dst_1.9x2.5_c090203.nc'); names(nc$var)
lon.Z2 = ncvar_get(nc, 'lon'); lat.Z2 = ncvar_get(nc, 'lat')
lon.Z2 = lon.Z2[c(73:144,1:72)]; lon.Z2[which(lon.Z2>=180)] = lon.Z2[which(lon.Z2>=180)] - 360
mbl_bsn_fct_Z2 = ncvar_get(nc,'mbl_bsn_fct_geo')[c(73:144,1:72),]
plot.field.log(mbl_bsn_fct_Z2, lon.Z2, lat.Z2, legend.mar=5)
mbl_bsn_fct_Z.05x06 = sp.smooth(mbl_bsn_fct_Z2, lon.Z2, lat.Z2, lon.M2, lat.M2)
plot.field(mbl_bsn_fct_Z.05x06, lon.M2, lat.M2)
plot.field.log(mbl_bsn_fct_Z.05x06, lon.M2, lat.M2, legend.mar=5)
save('mbl_bsn_fct_Z.05x06', 'lon.M2', 'lat.M2', file='zender_dust_source_05x0625.RData')

#############################################################
# intermittency and boundary-layer meteorology
nc = nc_open('/Volumes/SEAGATE/MERRA2/upscale_2006/tavg1_2d_lnd_Nx/MERRA2_300.tavg1_2d_lnd_Nx.20060101.SUB.nc')
SHLAND = ncvar_get(nc, 'SHLAND')[,,n] # LAI (W / m2)
#plot.field(SHLAND[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5)
#plot.field(SHLAND, lon.M2, lat.M2, legend.mar=5, type='sign', col=TEMP_DIFF_65, zlim=c(-500,500))
nc = nc_open('/Volumes/SEAGATE/MERRA2/upscale_2006/tavg1_2d_slv_Nx/MERRA2_300.tavg1_2d_slv_Nx.20060101.SUB.nc')
T10M = ncvar_get(nc, 'T10M')[,,n]
#plot.field(T10M, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(220,320))
nc = nc_open('/Volumes/SEAGATE/MERRA2/upscale_2006/tavg1_2d_flx_Nx/MERRA2_300.tavg1_2d_flx_Nx.20060101.SUB.nc'); names(nc$var)
USTAR = ncvar_get(nc, 'USTAR')[,,n]
RHOA = ncvar_get(nc, 'RHOA')[,,n]
PBLH = ncvar_get(nc, 'PBLH')[,,n]
Z0M = ncvar_get(nc, 'Z0M')[,,n]
plot.field(PBLH, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(0,4000), col=WBGYR)
plot.field(PBLH[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5, type='def', zlim=c(0,4000), col=WBGYR)
plot.field(USTAR[iC05-288,90:170], lon.M2[iC05-288], lat.M2[90:170], legend.mar=5, type='def', zlim=c(0,1))
plot.field.log(Z0M[iC05-288,90:170], lon.M2[iC05-288], lat.M2[90:170], legend.mar=5, type='def', zlim=c(-4,1))
plot.field(USTAR, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(0,1))

# Monin-Obukhov length (Donald Golder, 1972; Khaled Essa, 1999)
c_p = 1005   # J/kg/K  specific heat of air at constant pressure
k = 0.4      # von Karmen constant
g = 9.81     # m / s^2, gravity 
Obu_L = - RHOA * c_p * T10M * USTAR^3 / (k * g * SHLAND)
plot.field(-Obu_L[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5, type='def', col=WBGYR, zlim=c(0,400))
plot.field(Obu_L[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5, type='sign', col=TEMP_DIFF_65, zlim=c(-4000,4000))
plot.field(Obu_L, lon.M2, lat.M2, legend.mar=5, type='sign', col=TEMP_DIFF_65, zlim=c(-400,400))

plot.field(-0.5*(PBLH/Obu_L)[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5, type='def', col=WBGYR, zlim=c(0,400))
plot.field(12-0.5*(PBLH/Obu_L)[indi,indj], lon.M2[indi], lat.M2[indj], legend.mar=5, type='def', col=WBGYR, zlim=c(10,20))
plot.field((12-0.5*(PBLH/Obu_L)[indi,indj])^0.333, lon.M2[indi], lat.M2[indj], legend.mar=5, type='def', col=WBGYR, zlim=c(2.3,8))
12^0.333
plot.field(-0.5*(PBLH/Obu_L), lon.M2, lat.M2, legend.mar=5, type='sign', zlim=c(-40,40))

# construct intermittency variables
# first choose a friction velocity. Aeolian/Aerodynamic? Neutral/Non-neutral?
# USTAR.aeo.neu, USTAR.aeo.nonneu, USTAR.aero.neu, USTAR (aero+nonneu, default MERRA2 output)
wnd_frc_slt = USTAR* F_eff.LC0   #  default MERRA2 (aero+nonneu)output translated to 0.1 m saltation height
wnd_frc_slt = USTAR.aeo.neu*F_eff.LC0   #  default MERRA2 (aero+nonneu)output translated to 0.1 m saltation height
wnd_frc_slt = USTAR.aeo.nonneu*F_eff.LC0   #  default MERRA2 (aero+nonneu)output translated to 0.1 m saltation height
u_mean_slt = (wnd_frc_slt/k) * log(0.1 / 1e-4)  # We used z0a = 1e-4 m in CESM and assumed. Should we change to Pr05 here?
# sd of lowpass-filtered wind speed
#u_sd_slt = wnd_frc_slt * (12)^0.333   # note that there is no buoyancy term here; need to add in later when Obukhov length is obtained (1 Oct 2021)
u_sd_slt = wnd_frc_slt * (12 - 0.5*PBLH/Obu_L)^0.333   # buoyancy is added (4 Oct 2021)
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

plot.field(u_mean_slt, lon.M2, lat.M2)
plot.field(u_sd_slt, lon.M2, lat.M2)
plot.field(u_sd_slt[indi,indj], lon.M2[indi], lat.M2[indj], col=WBGYR, type='def', zlim=c(0,1))
plot.field(u_sd_slt_b[indi,indj], lon.M2[indi], lat.M2[indj], col=WBGYR, type='def', zlim=c(0,1))
plot.field(u_fld_thr, lon.M2, lat.M2)
plot.field(u_impct_thr, lon.M2, lat.M2)
plot.field(intrmtncy_fct, lon.M2, lat.M2, col=WBGYR)
plot.field.log(intrmtncy_fct, lon.M2, lat.M2, col=WBGYR, legend.mar=5, type='def', zlim=c(-7,0))
plot.field(intrmtncy_fct[indi,indj], lon.M2[indi], lat.M2[indj], col=WBGYR, legend.mar=5)
plot.field((eta_b-intrmtncy_fct)[indi,indj], lon.M2[indi], lat.M2[indj], type='def', zlim=c(0,0.2), legend.mar=5)
# what else is needed
# liquid fraction?

###################
# a few more essential variables for dust emission equations
# standardized fluid threshold (m/s)
rho_a0 = 1.225; # kg m-3, standard atmospheric density
ustar_st0 = 0.16;  # m s-1, minimum value of standardized fluid threshold
C_d0 = 4.4e-5; C_e = 2.0; C_alpha = 2.7 # dust emission coefficient parameters
ustar_st = ustar_ft.wet * sqrt(RHOA/rho_a0)
# dust emission coefficient (a measure of soil erodibility)
C_d = C_d0 * exp(-C_e * (ustar_st-ustar_st0)/ustar_st0);
#A = apply(C_d, c(1,2), mean, na.rm=T)
plot.field.log(C_d/C_d0, lon.M2, lat.M2, type='def', zlim=c(-5,0), legend.mar=5)

###################
# compute 3-D dust emissions
frag_expt_lim = 3    # fragmentation exponent limit (CESM simulations can blow up when 
# first get fragmentation exponent and limit to frag_expt_lim and below
frag_expt = (C_alpha*(ustar_st-ustar_st0)/ustar_st0)
frag_expt[which(frag_expt>frag_expt_lim)] = frag_expt_lim
A = apply(frag_expt, c(1,2), mean, na.rm=T)
plot.field(A, lon.M2, lat.M2, legend.mar=5)

# compute dust emissions using C19
C_tune = 0.05; dt = 2
F_d = array(0, dim=dim(C_d))   # using dimension of C_d to define F_d dimension
ind = which(wnd_frc_slt > ustar_it)  # 3-D index for which wind speed > threshold
F_d[ind] = C_tune * C_d[ind] * f_bare[ind] * RHOA[ind] * (wnd_frc_slt[ind]^2-ustar_it[ind]^2)/ustar_it[ind] * (wnd_frc_slt[ind]/ustar_it[ind])^(frag_expt[ind]) * intrmtncy_fct[ind]             # calculate F_d for indices with wind speed > threshold
F_d = sweep(F_d, c(1,2), mss_frc_cly_vld.FAO, '*')  # multiply 3-D fields by 2-D fields
F_d_sum = apply(F_d, c(1,2), sum, na.rm=T) * 3600 #* dt * 365 * 12
plot.field.log(F_d_sum, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(-7,-3))

# default K14
F_d = array(0, dim=dim(C_d))   # using dimension of C_d to define F_d dimension
ind = which(wnd_frc_slt > ustar_ft.wet)  # 3-D index for which wind speed > threshold
F_d[ind] = C_tune * C_d[ind] * f_bare[ind] * RHOA[ind] * (wnd_frc_slt[ind]^2-ustar_ft.wet[ind]^2)/ustar_ft.wet[ind] * (wnd_frc_slt[ind]/ustar_ft.wet[ind])^(frag_expt[ind])      # calculate F_d for indices with wind speed > threshold
F_d = sweep(F_d, c(1,2), mss_frc_cly_vld.FAO, '*')  # multiply 3-D fields by 2-D fields
F_d_sum = apply(F_d, c(1,2), sum, na.rm=T) * 3600 #* dt * 365 * 12
plot.field.log(F_d_sum, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(-7,-3))


#############################################################
#############################################################

# dust emission code trial: 3 days of MERRA2 data for testing
# extract N Africa for calculation only
# allow 1-hourly, 2-hourly data or other time resolution simulations (can skip some timestep and multiply dt at the end if one think e.g. 2-hourly sampling is good enough)
# allow showing actual dust emissions and annually scaled emissions (Tg / yr)
# set up spatial and temporal domain

# 30 Sep 2021
# remaining issues:
# 1. Confirm does CESM using porosity or effective porosity (discounting ice portion) to calculate gravimetric water content? We use MERRA2 here which uses porosity
# 2. CESM uses liquid water fractions for dust emissions, but here we haven't added this in for MERRA2 dust emissions
# 3. For intermittency, we assumed z0a = 1e-4 m globally. Change to hybrid Z0a? Are we going to do the same for CESM?
# 4. For intermittency, we used log law to translate friction velocity to 0. m saltation level. Do we need to add in non-neutral term?

##################
# user input here
# specialize temporal domain as YYYYMMDD
startdate = 20060101   # YYYYMMDD
enddate = 20060103   # YYYYMMDD
# specialize temporal resolution (dt = 1 as hourly, dt = 3 as 3-hourly etc.)
dt = 2   # try 2-hourly
# specialize spatial domain as indices
indi = 251:390; indj = 171:260     # lon and lat indices for N Africa
# choose clay dataset: 'FAO' or 'SG' (SoilGrids)
clay_dataset = 'FAO'

##################
# calculate model parameters based on user inputs
date_vec = make.date.vec(startdate, enddate)  # make date vec
time = 1:(length(date_vec)*24/dt)  # a vector of timesteps
indt = seq(dt, 24, by=dt)   # hour of day
# get
nc = nc_open('/Volumes/GoogleDrive/My Drive/MERRA2/upscale_2006/MERRA2_300.tavg1_2d_lnd_Nx.20060101.SUB.nc')
lon.M2 = ncvar_get(nc, 'lon')[indi]; lat.M2 = ncvar_get(nc, 'lat')[indj]
nc_close(nc)

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
frag_expt_lim = 7    # fragmentation exponent limit (CESM simulations can blow up when frag_expt is too large when using impact threshold for dust emission equation)

###################
# load required time-invariant datasets
.env = new.env()   # set a new environment for storing data temporarily
# load GLCNMO LULC data
load('GLCNMO_frc_area_06x0625.RData',envir=.env); ls(.env)
frc.r.M2 = .env$frc.r.M2[indi,indj]; frc.veg.M2 = .env$frc.veg.M2[indi,indj]
frc.mix.M2 = .env$frc.mix.M2[indi,indj]
A_r = frc.r.M2 / (frc.r.M2+frc.veg.M2+frc.mix.M2)
A_veg = frc.veg.M2 / (frc.r.M2+frc.veg.M2+frc.mix.M2)
A_mix = frc.mix.M2 / (frc.r.M2+frc.veg.M2+frc.mix.M2)
# load landfilt 
load('landfilt_05x0625.RData',envir=.env)
landfilt.M2 = .env$landfilt.M2[indi,indj]
# get Prigent 2005 roughness length and rock drag partition factor
load('Pr05_z0.min_05x0625.RData',envir=.env)
z0.r = .env$z0.05[indi,indj]/100   # in cm, divide by 100 to convert to m
# in this formulation, set globally constant z0s because soil median diameter ~ 130 um
f_eff.r = 1 - ( log(z0.r/z0s) / log(0.7*(X/z0s)^0.8) )  
f_eff.r = f_eff.r*landfilt.M2   # take away values in Caspian sea
# image(f_eff.r)
# get porosity  (issue: MERRA2 porosity may be different from CESM's effective porosity)
nc = nc_open('/Volumes/GoogleDrive/My Drive/MERRA2/upscale_2006/MERRA2_100.const_2d_lnd_Nx.00000000.nc4'); names(nc$var)
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
  }
  
  else if (clay_dataset == 'SG') {
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

###################
# scrap files for selected dates
files_flx = files_slv = files_lnd = NULL   # character strings storing the file names
for (dd in 1:length(date_vec)) {
  date_vec[dd]
  files_flx = c(files_flx, Sys.glob(paste('/Volumes/GoogleDrive/My Drive/MERRA2/upscale_2006/MERRA2*flx*', date_vec[dd], '*.nc', sep='')))  # for surface flux variables: surface USTAR, Z0M, and RHOA (air density)
  files_slv = c(files_slv, Sys.glob(paste('/Volumes/GoogleDrive/My Drive/MERRA2/upscale_2006/MERRA2*slv*', date_vec[dd], '*.nc', sep='')))  # for surface level variables: surface U10, V10, and DISPH (displacement height)
  files_lnd = c(files_lnd, Sys.glob(paste('/Volumes/GoogleDrive/My Drive/MERRA2/upscale_2006/MERRA2*lnd*', date_vec[dd], '*.nc', sep='')))  # for land variables: SFMC (Volumetric water content) and LAI
}

###################
# start the main loop
# for days in d
d = 1   # d-th day

###################
# get 3D data in time and space, subsetted as regional data with variable time resolutions
i1 = indi[1]; j1 = indj[1]; t1 = indt[1]
il = length(indi); jl = length(indj); tl = length(indt)

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
nc = nc_open(files_flx[d]); names(nc$var)#, start=c(11,26), count=c(5,5)
USTAR = ncvar_get(nc, 'USTAR', start=c(i1,j1,t1), count=c(il,jl,tl)) # m/s
Z0M = ncvar_get(nc, 'Z0M', start=c(i1,j1,t1), count=c(il,jl,tl))     # m
RHOA = ncvar_get(nc, 'RHOA', start=c(i1,j1,t1), count=c(il,jl,tl))     # kg/m3
PBLH = ncvar_get(nc, 'PBLH', start=c(i1,j1,t1), count=c(il,jl,tl))     # kg/m3
nc = nc_open(files_slv[d]); names(nc$var)
DISPH = ncvar_get(nc, 'DISPH', start=c(i1,j1,t1), count=c(il,jl,tl))   # m
U10M = ncvar_get(nc, 'U10M', start=c(i1,j1,t1), count=c(il,jl,tl))     # m/s
V10M = ncvar_get(nc, 'V10M', start=c(i1,j1,t1), count=c(il,jl,tl))     # m/s
T10M = ncvar_get(nc, 'T10M', start=c(i1,j1,t1), count=c(il,jl,tl))     # K
U10 = sqrt(U10M^2 + V10M^2)  # m/s. Define U10 as the total wind speed of U10M and V10M. Don't confuse U10 with U10M. 
#proc.time() - ptm

rm('DISPH', 'U10M', 'V10M') # remove data to save space
nc = nc_open(files_lnd[d]); names(nc$var)
SFMC = ncvar_get(nc, 'SFMC', start=c(i1,j1,t1), count=c(il,jl,tl))
LAI = ncvar_get(nc, 'LAI', start=c(i1,j1,t1), count=c(il,jl,tl))
SHLAND = ncvar_get(nc, 'SHLAND', start=c(i1,j1,t1), count=c(il,jl,tl))
#plot.field(DISPH[,,3], lon.M2, lat.M2)
#plot.field(USTAR[,,3], lon.M2, lat.M2)
#plot.field(Z0M[,,3], lon.M2, lat.M2)
#plot.map.log((Z0M[,,3]),lon.M2, lat.M2,type='def',zlim=c(-4,0), legend.mar=5)

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
emisfilt = f_bare/f_bare    # this filter is 3D and changes with time; or one can make a annual, 2D one. 
#image(emisfilt)

# Okin's vegetation drag partition
# Using Caroline Pierre (2014)'s formulation of Okin's scheme
K = pi/2 * (1/LAI - 1); K[which(K<0)] = 0
SSR = (K + U_0 * c) / (K + c)  # Pierre-Okin shear stress ratio (SSR)
#plot.field(SSR[,,3], lon.M2, lat.M2, type='def', zlim=c(0.6,0.9))

# combine rock and vegetation drag partition factors
# then combine drag partition factors
F_eff.LC1 = array(NaN, dim=c(length(indi), length(indj), length(indt)))
#F_eff.LC2 = array(NaN, dim=c(length(indi), length(indj), length(indt)))
for (h in 1:length(indt)) {
  #F_eff.LC1[,,h] = (frc.r.M2*f_eff.r^3 + frc.veg.M2*SSR[,,h]^3 + frc.mix.M2*(SSR[,,h]*f_eff.r)^3)^(1/3)
  F_eff.LC1[,,h] = (A_r*f_eff.r^3 + A_veg*SSR[,,h]^3 + A_mix*(SSR[,,h]*f_eff.r)^3)^(1/3)
  #F_eff.LC2[,,h] = (frc.r.M2*f_eff.r^3 + (frc.veg.M2+frc.mix.M2)*SSR[,,h]^3)^(1/3)  # newly define LC2, removed mixed regime
}

#plot.field(F_eff.LC1[,,6], lon.M2, lat.M2, type='def', zlim=c(0,1), col=WBGYR)
#plot.field(F_eff.LC2[,,6], lon.M2, lat.M2, type='def', zlim=c(0,1), col=WBGYR)

# then get hybrid aeolian roughness length (Z0.LC1) from hybrid F_eff, using Marticorena's equation
#f_eff.r = 1 - ( log(z0.r/z0s) / log(0.7*(X/z0s)^0.8) )  
Z0.LC1 = z0s * exp( (1 - F_eff.LC1) * log(0.7*(X/z0s)^0.8) )   # in m
#plot.field.log(Z0.LC1[,,6], lon.M2, lat.M2, type='def', zlim=c(-5,0), col=tim.colors(32)[32:1], legend.mar=5)
#plot.field.log(Z0M[,,6], lon.M2, lat.M2, type='def', zlim=c(-5,0), col=tim.colors(32)[32:1], legend.mar=5)
#plot.field.log(z0.r*landfilt.M2, lon.M2, lat.M2, type='def', zlim=c(-5,0), col=tim.colors(32)[32:1], legend.mar=5)   # also in m

###################
# reconstruct friction velocity with different meanings and choose which to use

# test if U10, USTAR, and Z0M are consistent with each other
KSI_m = k*U10/USTAR - log(10 / Z0M) # first determine MOST instability (KSI)
for (aa in 1) {
  if (wind.scale=='default') {print('default MERRA2 ustar is used')}
  else if (wind.scale=='aeo_nonneu') {USTAR.aeo.nonneu = k * U10 / (log(10 / Z0.LC1)+KSI_m); print('aeolian+non-neutral ustar is used')}  # reconstruct USTAR from U10 and Z0M
  #USTAR.aero.neu = k * U10 / log(10 / Z0M)  # reconstruct USTAR from U10 and Z0M; aero means aerodynamic, neu means neutral (log wind profile) 
  # calculate aeolian ustar assuming neutral condition
  else if (wind.scale=='aeo_neu') {USTAR.aeo.neu = k * U10 / log(10 / Z0.LC1); print('aeolian+neutral ustar is used')}  # reconstruct USTAR from U10 and Z0M
}
#plot.field(USTAR[,,6]*landfilt.M2, lon.M2, lat.M2, type='def', zlim=c(0,0.8), col=WBGYR)
#plot.field(USTAR.aeo.nonneu[,,6], lon.M2, lat.M2, type='def', zlim=c(0,0.8), col=WBGYR)
#plot.field(USTAR.aeo.nonneu[,,6]/USTAR.aeo.neu[,,6], lon.M2, lat.M2, type='def', zlim=c(0,2), col=TEMP_DIFF_65)
#plot.field(USTAR.aeo.neu[,,6], lon.M2, lat.M2, type='def', zlim=c(0,0.8), col=WBGYR)
A = apply(USTAR, c(1,2), mean, na.rm=T)
plot.field(A*landfilt.M2, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(0,0.7))
A = apply(USTAR.aeo.nonneu, c(1,2), mean, na.rm=T)
plot.field(A, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(0,0.7))
A = apply(USTAR.aeo.neu, c(1,2), mean, na.rm=T)
plot.field(A, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(0,0.7))

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
A = apply(USTAR - ustar_it, c(1,2), mean, na.rm=T)
plot.field(A*landfilt.M2, lon.M2, lat.M2, legend.mar=5, type='sign', zlim=c(-0.5,0.5), col=TEMP_DIFF_65)
A = apply(USTAR.aeo.nonneu - ustar_it, c(1,2), mean, na.rm=T)
plot.field(A, lon.M2, lat.M2, legend.mar=5, type='sign', zlim=c(-0.5,0.5), col=TEMP_DIFF_65)

###################
# construct intermittency

# first get Monin-Obukhov length in m (Donald Golder, 1972; Khaled Essa, 1999)
Obu_L = - RHOA * c_p * T10M * USTAR^3 / (k * g * SHLAND)

# then choose a friction velocity. Aeolian/Aerodynamic? Neutral/Non-neutral?
# USTAR.aeo.neu, USTAR.aeo.nonneu, USTAR.aero.neu, USTAR (aero+nonneu, default MERRA2 output)
for (aa in 1) {
  if (wind.scale=='default') {wnd_frc_slt = USTAR*F_eff.LC1; print('default MERRA2 ustar is used for intermittency and emission')}   #  default MERRA2 (aero+nonneu)output
  else if (wind.scale=='aeo+neu') {wnd_frc_slt = USTAR.aeo.neu*F_eff.LC1; print('aeolian+neutral ustar is used for intermittency and emission')}   #  reconstructed aeo+neu USTAR
  else if (wind.scale=='aeo+nonneu') {wnd_frc_slt = USTAR.aeo.nonneu*F_eff.LC1; print('aeolian+neutral ustar is used for intermittency and emission')}   #  reconstructed aeo+nonneu USTAR
}
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

plot.field(u_mean_slt[,,6], lon.M2, lat.M2)
plot.field(u_sd_slt[,,6], lon.M2, lat.M2)
plot.field(u_fld_thr[,,6], lon.M2, lat.M2)
plot.field(u_impct_thr[,,6], lon.M2, lat.M2)
plot.field(intrmtncy_fct[,,6], lon.M2, lat.M2)
A = apply(intrmtncy_fct, c(1,2), mean, na.rm=T)
plot.field(A, lon.M2, lat.M2, col=WBGYR, legend.mar=5, type='def', zlim=c(0,1))

# plot wind time series
#ii = 23:25; jj = 47:50  # el djouf
ii = 24; jj = 48
tt = 1:12*2
par(mai=c(0.6,0.6,0.1,0.1))
plot(NA,NA, xlim=c(1,24),ylim=c(0,10), xlab='hour UTC+0', ylab='wind speed (m/s)')
lines(tt, u_fld_thr[ii,jj,], 'l', col='black', lwd=2)
lines(tt, u_impct_thr[ii,jj,], 'l', col='black', lwd=2)
lines(tt, u_mean_slt[ii,jj,], 'l', col='blue', lwd=2)
lines(tt, (u_mean_slt+2*u_sd_slt)[ii,jj,], 'l', col='red', lwd=2, lty=2)
lines(tt, (u_mean_slt-2*u_sd_slt)[ii,jj,], 'l', col='red', lwd=2, lty=2)


###################
# a few more essential variables for dust emission equations
# standardized fluid threshold (m/s)
ustar_st = ustar_ft.wet * sqrt(RHOA/rho_a0)
# dust emission coefficient (a measure of soil erodibility)
C_d = C_d0 * exp(-C_e * (ustar_st-ustar_st0)/ustar_st0);
A = apply(C_d, c(1,2), mean, na.rm=T)
plot.field.log(A/C_d0, lon.M2, lat.M2, type='def', zlim=c(-5,0), legend.mar=5)


###################
# compute 3-D dust emissions

# first get fragmentation exponent and limit to frag_expt_lim and below
frag_expt = (C_alpha*(ustar_st-ustar_st0)/ustar_st0)
frag_expt[which(frag_expt>frag_expt_lim)] = frag_expt_lim
A = apply(frag_expt, c(1,2), mean, na.rm=T)
plot.field(A, lon.M2, lat.M2, legend.mar=5)

# compute dust emissions
F_d = array(0, dim=dim(C_d))   # using dimension of C_d to define F_d dimension
ind = which(wnd_frc_slt > ustar_it)  # 3-D index for which wind speed > threshold
F_d[ind] = C_tune * C_d[ind] * f_bare[ind] * RHOA[ind] * (wnd_frc_slt[ind]^2-ustar_it[ind]^2)/ustar_it[ind] * (wnd_frc_slt[ind]/ustar_it[ind])^(frag_expt[ind]) * intrmtncy_fct[ind]             # calculate F_d for indices with wind speed > threshold
F_d = sweep(F_d, c(1,2), mss_frc_cly_vld, '*')  # multiply 3-D fields by 2-D fields
F_d_sum = apply(F_d, c(1,2), sum, na.rm=T) * 3600 * dt * 365
plot.field.log(F_d_sum, lon.M2, lat.M2, legend.mar=5, col=WBGYR)

filename = paste('dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')
save('F_d', 'lon.M2', 'lat.M2', file=)



