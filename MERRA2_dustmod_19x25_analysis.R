#############################################################
#############################################################

# 12 Oct 2021
# dust emission code analysis
# functions for plot and analysis

# set working directory
setwd('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR')
# load libraries 
library(ncdf4); library(fields); library(maps); library(Metrics); library(pracma); library(abind)
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
indi = 1:144; indj = 18:93 
# choose clay dataset: 'FAO' or 'SG' (SoilGrids)
clay_dataset = 'FAO'
# choose the wind scale for dust emission, including default (aerodynamic + nonneutral) ustar (default), aeolian + nonneutral star (aeo+nonneu), and aeolian + neutral ustar (aeo+neu)
wind.scale='default'
# choose land cover averaging method: LC0 (rock and veg) or LC1 (rock, veg and mix)
LC = 'LC0'
# choose the averaging method of the coarse input variables: 'mean' or 'aw_mean' (area-weighted mean)
avg.method = 'aw_mean'
# choose whether default Kok et al. (2014) scheme (K14) or intermittency (K14+Comola et al., 2019) is used. 'FALSE' is K14, 'TRUE' is K14+C19 intermittency scheme.
use.intermittency = TRUE
# fragmentation exponent limit (CESM simulations can blow up when frag_expt is too large when using impact threshold for dust emission equation)
frag_expt_lim = 3   
# dust production model (DPM): Kok et al. (2014) K14 or Zender et al. (2003) Z03
DPM = 'Z03'
if (DPM=='Z03') use.intermittency = 'FALSE' # by default Z03 has no intermittency eff
# if source function is needed (for Z03), choose 'Zender' or 'Ginoux'. For K14, type 'FALSE'
src_fn = 'Ginoux'
##################
# calculate model parameters based on user inputs
date_vec = make.date.vec(startdate, enddate)  # make date vec
time = 1:(length(date_vec)*24/dt)  # a vector of timesteps
indt = seq(dt, 24, by=dt)   # hour of day
# get CESM lon lat
files.clm = Sys.glob('/Volumes/GoogleDrive/My Drive/CESM_tempfiles/renewmod5_20042005/new_intermittency/*.clm*.nc'); nc = nc_open(files.clm[1])
LON.CESM = seq(-180,177.5,by=2.5)[indi]; LAT.CESM = ncvar_get(nc, 'lat')[indj]
indlon = c(73:144,1:72); LON.r = LON.CESM[indlon]; LON.r[141:144]=LON.r[141:144]+360
nc_close(nc)

##################
# get a 2-D area map for MERRA2 resolution
get.area = function(lon, lat, unit='km2') {
  # allow two units: km2 (default) and m2.
  dlon = lon[2]-lon[1]; dlat = lat[2]-lat[1]
  lat.area = NULL
  for (j in 1:length(lat)) {
    lat1 = lat[j]-dlat/2; lat2 = lat[j]+dlat/2
    lon1 = lon[j]-dlon/2; lon2 = lon[j]+dlon/2
    lat.area = c(lat.area, area.latlon(lat1, lon1, lat2, lon2))
  }
  area = NULL; for (i in 1:length(lon)) {area = rbind(area, lat.area)}
  if (unit=='m2') area = area*1e6  # from km2 to m2
  return(area)
}

area.CESM = get.area(LON.CESM, LAT.CESM, 'm2')
plot.field(area.CESM, LON.CESM, LAT.CESM, legend.mar=5)

##################
# look at the sum of gridded dust emissions
.env = new.env()
v = 1:length(indt)
#files = Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=3/dust_emis_flx_19x25_*.RData', sep=''))
files.int = Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=TRUE/dust_emis_flx_19x25_*.RData', sep=''))
files.noint = Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=FALSE/dust_emis_flx_19x25_*.RData', sep=''))
F_d.noint = F_d.int = array(NaN, dim=c(length(LON.CESM), length(LAT.CESM), length(files)*length(indt)))
for (n in 1:length(files)) {
  load(files.int[n],envir=.env); F_d.int[,,v] = .env$F_d
  load(files.noint[n],envir=.env); F_d.noint[,,v] = .env$F_d 
  v = v+length(indt)
}
F_d.CESM.int = apply(F_d.int, c(1,2), sum, na.rm=T) * 3600 * dt 
sum(F_d.CESM.int*area.CESM, na.rm=TRUE)/1e9
F_d.CESM.noint = apply(F_d.noint, c(1,2), sum, na.rm=T) * 3600 * dt 
sum(F_d.CESM.noint*area.CESM, na.rm=TRUE)/1e9
plot.field.log(F_d.CESM.int, LON.CESM, LAT.CESM, legend.mar=5, col=WBGYR, type='def', zlim=c(-4,0))
plot.field.log(F_d.CESM.noint, LON.CESM, LAT.CESM, legend.mar=5, col=WBGYR, type='def', zlim=c(-4,0))
plot.field.log(F_d.CESM.int, LON.CESM, LAT.CESM, legend.mar=5, type='def', zlim=c(-4,0))
plot.field.log(F_d.CESM.noint, LON.CESM, LAT.CESM, legend.mar=5, type='def', zlim=c(-4,0))

# plot aggregated 0.5x0.625 dust emissions
files = Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe/dust_emis_flx_05x0625_*.RData', sep=''))
.env1 = new.env()
for (d in 1:length(files)) {
  print(d)
  load(files[d], envir=.env1)
  if (d==1) {lon.M2 = .env1$lon.M2; lat.M2 = .env1$lat.M2; F_d.M2 = .env1$F_d}
  else if (d >= 2) {F_d.M2 = F_d.M2 + .env1$F_d}
  #F_d = abind(F_d, (apply(.env1$F_d, c(1,2), sum, na.rm=T)*3600*dt), along=3)
  #F_d = apply(F_d, c(1,2), sum, na.rm=T)
}
dim(F_d.M2)
F_d.M2 = apply(F_d.M2, c(1,2), sum, na.rm=T)*3600*dt
plot.field.log(F_d.M2, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(-4,1))
plot.field(F_d.M2, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(0,3), col=WBGYR)
area.M2 = get.area(lon.M2, lat.M2, 'm2')#; rm(area.M2)
sum(F_d.M2*area.M2, na.rm=TRUE)/1e9
plot.field(area.M2, lon.M2, lat.M2, legend.mar=5)
F_d.M2.19x25 = sp.dissolve(F_d.M2, lon.M2, lat.M2, LON.CESM, LAT.CESM)
plot.field.log(F_d.M2.19x25, LON.CESM, LAT.CESM, legend.mar=5, type='def', zlim=c(-4,1))
sum(F_d.M2.19x25*area.CESM, na.rm=TRUE)/1e9

R = F_d.M2.19x25/F_d.CESM; R[which(is.infinite(R) | R>(1e3) | R<(1e-3))]=NaN
D = F_d.M2.19x25-F_d.CESM
plot.field.log(R, LON.CESM, LAT.CESM, legend.mar=5, type='def', zlim=c(-3,3), col=TEMP_DIFF_65)
plot.field(D, LON.CESM, LAT.CESM, legend.mar=5, type='sign', zlim=c(-1,1), col=TEMP_DIFF_65)


###################
# do the same thing for Z03 scheme
# added 18 Nov 2021
.env = new.env()
v = 1:length(indt)
files = Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_',DPM,'_src_fn=',src_fn,'/dust_emis_flx_19x25_*.RData', sep=''))
F_d.Z03.G.CESM = array(NaN, dim=c(length(LON.CESM), length(LAT.CESM), length(files)*length(indt)))
for (n in 1:length(files)) {
  load(files[n],envir=.env); F_d.Z03.G.CESM[,,v] = .env$F_d
  v = v+length(indt)
}
F.Z03.G.CESM.avg = apply(F_d.Z03.G.CESM, c(1,2), sum, na.rm=T) * 3600 * dt 
sum(F.Z03.G.CESM.avg*area.CESM, na.rm=TRUE)/1e9
plot.field.log(F.Z03.G.CESM.avg, LON.CESM, LAT.CESM, legend.mar=5, col=WBGYR, type='def', zlim=c(-4,0))
plot.field.log(F.Z03.G.CESM.avg, LON.CESM, LAT.CESM, legend.mar=5, type='def', zlim=c(-4,0))

# plot aggregated 0.5x0.625 dust emissions
files = Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe/dust_emis_flx_05x0625_*.RData', sep=''))
.env1 = new.env()
for (d in 1:length(files)) {
  print(d)
  load(files[d], envir=.env1)
  if (d==1) {lon.M2 = .env1$lon.M2; lat.M2 = .env1$lat.M2; F_d.M2 = .env1$F_d}
  else if (d >= 2) {F_d.M2 = F_d.M2 + .env1$F_d}
  #F_d = abind(F_d, (apply(.env1$F_d, c(1,2), sum, na.rm=T)*3600*dt), along=3)
  #F_d = apply(F_d, c(1,2), sum, na.rm=T)
}
dim(F_d.M2)
F_d.M2 = apply(F_d.M2, c(1,2), sum, na.rm=T)*3600*dt
plot.field.log(F_d.M2, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(-4,1))
plot.field(F_d.M2, lon.M2, lat.M2, legend.mar=5, type='def', zlim=c(0,3), col=WBGYR)
area.M2 = get.area(lon.M2, lat.M2, 'm2')#; rm(area.M2)
sum(F_d.M2*area.M2, na.rm=TRUE)/1e9
plot.field(area.M2, lon.M2, lat.M2, legend.mar=5)
F_d.M2.19x25 = sp.dissolve(F_d.M2, lon.M2, lat.M2, LON.CESM, LAT.CESM)
plot.field.log(F_d.M2.19x25, LON.CESM, LAT.CESM, legend.mar=5, type='def', zlim=c(-4,1))
sum(F_d.M2.19x25*area.CESM, na.rm=TRUE)/1e9

R = F_d.M2.19x25/F_d.CESM; R[which(is.infinite(R) | R>(1e3) | R<(1e-3))]=NaN
D = F_d.M2.19x25-F_d.CESM
plot.field.log(R, LON.CESM, LAT.CESM, legend.mar=5, type='def', zlim=c(-3,3), col=TEMP_DIFF_65)
plot.field(D, LON.CESM, LAT.CESM, legend.mar=5, type='sign', zlim=c(-1,1), col=TEMP_DIFF_65)


###################
# get and plot input meteorology and land-surface variables

# first scrap files 
files = NULL   # character strings storing the file names
for (dd in 1:length(date_vec)) {
  date_vec[dd]
  files = c(files, Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/REGRIDDED_MERRA2/',avg.method,'/merra2_regridded_met_19x25_', date_vec[dd], '.RData', sep='')))  
}

U_0 = 0.32  # SSR in the immediate lee of a plant for Okin's scheme, dimensionless
c = 4.8     # e-folding distance SSR recovery for Okin's scheme, dimensionless
SFMC = LAI = SSR = SSR2 = NULL#array(NaN, dim=c(length(LON.CESM), length(LAT.CESM), length(indt)))
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
  if (d==1) {
  #LAI   = .env$LAI.CESM[i1:il,j1:jl,]
  # Okin's vegetation drag partition
  # Using Caroline Pierre (2014)'s formulation of Okin's scheme
  #K = pi/2 * (1/LAI - 1); K[which(K<0)] = 0
  #SSR = (K + U_0 * c) / (K + c)  # Pierre-Okin shear stress ratio (SSR)
  SFMC  = .env$SFMC.CESM[i1:il,j1:jl,]
  } else {
  #LAI   = LAI + .env$LAI.CESM[i1:il,j1:jl,]
  # Okin's vegetation drag partition
  # Using Caroline Pierre (2014)'s formulation of Okin's scheme
  #K = pi/2 * (1/.env$LAI.CESM[i1:il,j1:jl,] - 1); K[which(K<0)] = 0
  #SSR =  SSR + (K + U_0 * c) / (K + c)  # Pierre-Okin shear stress ratio (SSR)
  SFMC  = SFMC + .env$SFMC.CESM[i1:il,j1:jl,]
  }
}

LAI.avg = apply(LAI/365, c(1,2), mean, na.rm=T)
plot.field(LAI.avg, LON.CESM, LAT.CESM, type='def', zlim=c(0,1))
i = 60:100; j = 30:50
plot.field(LAI.avg[i,j], LON.CESM[i], LAT.CESM[j], type='def', zlim=c(0,1))
SSR.avg = apply(SSR/365, c(1,2), mean, na.rm=T)
SSR.avg[which(SSR.avg>0.97)]=NaN
plot.field(SSR.avg, LON.CESM, LAT.CESM, type='def', zlim=c(0,1))

max(SSR.avg, na.rm=T)

