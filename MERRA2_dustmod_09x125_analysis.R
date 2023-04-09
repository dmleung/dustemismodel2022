#############################################################
#############################################################

# 28 Oct 2021
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
use.intermittency = TRUE
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

area.C1 = get.area(LON.C1, LAT.C1, 'm2')
plot.field(area.C1, LON.C1, LAT.C1, legend.mar=5)

##################
# look at the sum of gridded dust emissions
.env = new.env()
v = 1:length(indt)
#files = Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=3/dust_emis_flx_1x1_*.RData', sep=''))
files.int = Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=TRUE/dust_emis_flx_1x1_*.RData', sep=''))
files.noint = Sys.glob(paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=FALSE/dust_emis_flx_1x1_*.RData', sep=''))
F_d.noint = F_d.int = array(NaN, dim=c(length(LON.C1), length(LAT.C1), length(files.noint)*length(indt)))
for (n in 1:length(files.noint)) {
  print(n)
  load(files.int[n],envir=.env); F_d.int[,,v] = .env$F_d
  load(files.noint[n],envir=.env); F_d.noint[,,v] = .env$F_d 
  v = v+length(indt)
}
F_d.C1.int = apply(F_d.int, c(1,2), sum, na.rm=T) * 3600 * dt 
sum(F_d.C1.int*area.C1, na.rm=TRUE)/1e9
F_d.C1.noint = apply(F_d.noint, c(1,2), sum, na.rm=T) * 3600 * dt 
sum(F_d.C1.noint*area.C1, na.rm=TRUE)/1e9
plot.field.log(F_d.C1.int, LON.C1, LAT.C1, legend.mar=5, col=WBGYR, type='def', zlim=c(-4,0))
plot.field.log(F_d.C1.noint, LON.C1, LAT.C1, legend.mar=5, col=WBGYR, type='def', zlim=c(-4,0))
plot.field.log(F_d.C1.int, LON.C1, LAT.C1, legend.mar=5, type='def', zlim=c(-4,0))
plot.field.log(F_d.C1.noint, LON.C1, LAT.C1, legend.mar=5, type='def', zlim=c(-4,0))

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
F_d.M2.1x1 = sp.dissolve(F_d.M2, lon.M2, lat.M2, LON.C1, LAT.C1)
plot.field.log(F_d.M2.1x1, LON.C1, LAT.C1, legend.mar=5, type='def', zlim=c(-4,1))
sum(F_d.M2.1x1*area.C1, na.rm=TRUE)/1e9

R = F_d.M2.1x1/F_d.C1; R[which(is.infinite(R) | R>(1e3) | R<(1e-3))]=NaN
D = F_d.M2.1x1-F_d.C1
plot.field.log(R, LON.C1, LAT.C1, legend.mar=5, type='def', zlim=c(-3,3), col=TEMP_DIFF_65)
plot.field(D, LON.C1, LAT.C1, legend.mar=5, type='sign', zlim=c(-1,1), col=TEMP_DIFF_65)
