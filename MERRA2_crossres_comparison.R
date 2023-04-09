#############################################################
#############################################################

# 18 Oct 2021
# cross-resolution dust emissions comparison
# functions for plot and analysis

# set working directory
setwd('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR')
# load libraries 
library(ncdf4); library(fields); library(maps); library(Metrics); library(pracma); library(abind); library(RColorBrewer)
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
indi.M2 = 1:576; indj.M2 = 65:350 
indi.C1 = 1:288; indj.C1 = 36:186 
indi.CESM = 1:144; indj.CESM = 18:93
indi.GC = 1:72; indj.GC = 9:44 
# choose clay dataset: 'FAO' or 'SG' (SoilGrids)
clay_dataset = 'FAO'
# choose the wind scale for dust emission, including default (aerodynamic + nonneutral) ustar (default), aeolian + nonneutral star (aeo+nonneu), and aeolian + neutral ustar (aeo+neu)
wind.scale='default'
# choose land cover averaging method: LC0 (rock and veg) or LC1 (rock, veg and mix)
LC = 'LC0'
# choose whether default Kok et al. (2014) scheme (K14) or intermittency (K14+Comola et al., 2019) is used. 'FALSE' is K14, 'TRUE' is K14+C19 intermittency scheme.
use.intermittency = TRUE
# fragmentation exponent limit (CESM simulations can blow up when frag_expt is too large when using impact threshold for dust emission equation)
frag_expt_lim = 3    
##################
# calculate model parameters based on user inputs
date_vec = make.date.vec(startdate, enddate)  # make date vec
time = 1:(length(date_vec)*24/dt)  # a vector of timesteps
indt = seq(dt, 24, by=dt)   # hour of day
# get MERRA2 lon lat
nc = nc_open('/Volumes/GoogleDrive/My Drive/MERRA2/upscale_2006/MERRA2_300.tavg1_2d_lnd_Nx.20060101.SUB.nc')
lon.M2 = ncvar_get(nc, 'lon')[indi.M2]; lat.M2 = ncvar_get(nc, 'lat')[indj.M2]
nc_close(nc)
i.M2 = c(20:576,1:19)
lon.M2.r = lon.M2[i.M2]; lon.M2.r[558:576]=lon.M2.r[558:576]+360
# load 0.9x1.25 CESM lat lon
nc = nc_open('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR/griddata_0.9x1.25_070212.nc')
LON.C1 = seq(-179.375,179.375,by=1.25)[indi.C1]; LAT.C1 = ncvar_get(nc, 'LATS')[1,2:192][indj.C1]
i.C1 = lon.ind = c(10:288,1:9); LON.C1.r = LON.C1[i.C1]; LON.C1.r[which(i.C1<=9)] = LON.C1.r[which(i.C1<=9)] + 360
# get CESM lon lat
files.clm = Sys.glob('/Volumes/GoogleDrive/My Drive/CESM_tempfiles/renewmod5_20042005/new_intermittency/*.clm*.nc'); nc = nc_open(files.clm[1])
LON.CESM = seq(-180,177.5,by=2.5)[indi.CESM]; LAT.CESM = ncvar_get(nc, 'lat')[indj.CESM]
indlon = c(6:144,1:5); LON.r = LON.CESM[indlon]; LON.r[140:144]=LON.r[140:144]+360
nc_close(nc)
# load 4x5 GC (GMAO) lat lon
LON.GC = seq(-180,175,by=5)[indi.GC]; LAT.GC = seq(-89,89,by=4)[indj.GC]
i.GC = c(4:72,1:3); LON.GC.r = LON.GC[i.GC]; LON.GC.r[70:72]=LON.GC.r[70:72]+360
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

area.M2 = get.area(lon.M2, lat.M2, 'm2')
#plot.field(area.M2, lon.M2, lat.M2, legend.mar=5)
area.C1 = get.area(LON.C1, LAT.C1, 'm2')
#plot.field(area.C1, LON.C1, LAT.C1, legend.mar=5)
area.CESM = get.area(LON.CESM, LAT.CESM, 'm2')
#plot.field(area.CESM, LON.CESM, LAT.CESM, legend.mar=5)
area.GC = get.area(LON.GC, LAT.GC, 'm2')
#plot.field(area.GC, LON.GC, LAT.GC, legend.mar=5)

# function for getting annual total dust emissions
get.tot.dustemis = function(date_vec, dt, res=NULL, use.intermittency=NULL) {
  # get annual dust emission total per grid (kg m-2 yr-1) 
  .env1 = new.env()
  for (d in 1:length(date_vec)) {
    print(date_vec[d])
    filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/dust_emis_flx_',res,'_', date_vec[d], '.RData', sep='')
    #filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006_new/', directory, '/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')
    load(filename, envir=.env1)
    if (d==1) {F_d=.env1$F_d; F_d[which(is.na(F_d))]=0}
    else if (d >= 2) {F_d1=.env1$F_d; F_d1[which(is.na(F_d1))]=0; F_d=F_d+F_d1}
    #F_d = abind(F_d, (apply(.env1$F_d, c(1,2), sum, na.rm=T)*3600*dt), along=3)
    #F_d = apply(F_d, c(1,2), sum, na.rm=T)
  }
  dim(F_d)
  F_d = apply(F_d, c(1,2), sum, na.rm=T)*3600*dt
  if (res=='05x0625') print(paste('total global dust emission is ',sum(F_d*area.M2, na.rm=TRUE)/1e9, ' Tg/yr',sep=''))
  else if (res=='19x25') print(paste('total global dust emission is ',sum(F_d*area.CESM, na.rm=TRUE)/1e9, ' Tg/yr',sep=''))
  return(list(F_d=F_d))
}


# extract annual dust emissions across different resolutions

#res = '05x0625' # First do MERRA2 resolution
# with intermittency
FF = get.tot.dustemis(date_vec, dt, res='05x0625', use.intermittency='TRUE'); F_d.M2.int = FF$F_d
sum(F_d.M2.int*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.M2.int[i.M2,]*4.5, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR); frame.region(frame='source', lwd=1.5)
F_d.M2.int.reg = extract.value.region(F_d.M2.int*area.M2/1e9, lon.M2, lat.M2, metric='sum', frame='source'); sum(F_d.M2.int*area.M2/1e9); sum(F_d.M2.int.reg)
# without intermittency
FF = get.tot.dustemis(date_vec, dt, res='05x0625', use.intermittency=FALSE); F_d.M2.noint = FF$F_d
plot.field.log(F_d.M2.noint, lon.M2, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR); frame.region(frame='source', lwd=1.5)
F_d.M2.noint.reg = extract.value.region(F_d.M2.noint*area.M2/1e9, lon.M2, lat.M2, metric='sum', frame='source'); sum(F_d.M2.noint*area.M2/1e9); sum(F_d.M2.noint.reg)

#res = '1x1' # CESM 1 deg resolution, actually 0.9x1.25 here
# with intermittency
FF = get.tot.dustemis(date_vec, dt, res='1x1', use.intermittency=TRUE); F_d.C1.int = FF$F_d
plot.field.log(F_d.C1.int[i.C1,]*4.5, LON.C1.r, LAT.C1, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR); frame.region(frame='source', lwd=1.5)
F_d.C1.int.reg = extract.value.region(F_d.C1.int*area.C1/1e9, LON.C1, LAT.C1, metric='sum', frame='source'); sum(F_d.C1.int*area.C1/1e9); sum(F_d.C1.int.reg)
# without intermittency
FF = get.tot.dustemis(date_vec, dt, res='1x1', use.intermittency=FALSE); F_d.C1.noint = FF$F_d
plot.field.log(F_d.C1.noint, LON.C1, LAT.C1, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR); frame.region(frame='source', lwd=1.5)
F_d.C1.noint.reg = extract.value.region(F_d.C1.noint*area.C1/1e9, LON.C1, LAT.C1, metric='sum', frame='source'); sum(F_d.C1.noint*area.C1/1e9); sum(F_d.C1.noint.reg)

#res = '19x25' # CESM 2 deg resolution
# with intermittency
FF = get.tot.dustemis(date_vec, dt, res='19x25', use.intermittency=TRUE); F_d.CESM.int = FF$F_d
plot.field.log(F_d.CESM.int[indlon,]*4.5, LON.r, LAT.CESM, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR); frame.region(frame='source', lwd=1.5)
F_d.CESM.int.reg = extract.value.region(F_d.CESM.int*area.CESM/1e9, LON.CESM, LAT.CESM, metric='sum', frame='source'); sum(F_d.CESM.int*area.CESM/1e9); sum(F_d.CESM.int.reg)
# without intermittency
FF = get.tot.dustemis(date_vec, dt, res='19x25', use.intermittency=FALSE); F_d.CESM.noint = FF$F_d
plot.field.log(F_d.CESM.noint, LON.CESM, LAT.CESM, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR); frame.region(frame='source', lwd=1.5)
F_d.CESM.noint.reg = extract.value.region(F_d.CESM.noint*area.CESM/1e9, LON.CESM, LAT.CESM, metric='sum', frame='source'); sum(F_d.CESM.noint*area.CESM/1e9); sum(F_d.CESM.noint.reg)

#res = '4x5' # GC (GMAO) 4x5 deg resolution
# with intermittency
FF = get.tot.dustemis(date_vec, dt, res='4x5', use.intermittency=TRUE); F_d.GC.int = FF$F_d
plot.field.log(F_d.GC.int*4.5, LON.GC, LAT.GC, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR); frame.region(frame='source', lwd=1.5)
F_d.GC.int.reg = extract.value.region(F_d.GC.int*area.GC/1e9, LON.GC, LAT.GC, metric='sum', frame='source'); sum(F_d.GC.int*area.GC/1e9); sum(F_d.GC.int.reg)
# without intermittency
FF = get.tot.dustemis(date_vec, dt, res='4x5', use.intermittency=FALSE); F_d.GC.noint = FF$F_d
plot.field.log(F_d.GC.noint, LON.GC, LAT.GC, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR); frame.region(frame='source', lwd=1.5)
F_d.GC.noint.reg = extract.value.region(F_d.GC.noint*area.GC/1e9, LON.GC, LAT.GC, metric='sum', frame='source'); sum(F_d.GC.noint*area.GC/1e9); sum(F_d.GC.noint.reg)

# Create a data frame
#AA = rep(0,time=9)
#MAT = rbind(F_d.M2.int.reg,F_d.C1.int.reg,F_d.CESM.int.reg,F_d.GC.int.reg,F_d.M2.noint.reg,F_d.C1.noint.reg,F_d.CESM.noint.reg,F_d.GC.noint.reg)
MAT = rbind(F_d.M2.int.reg,F_d.C1.int.reg,F_d.CESM.int.reg,F_d.GC.int.reg)*4.5
#col1 <- MAT[,1]; col2 <- MAT[,2]; col3 <- MAT[,3]; col4 <- MAT[,4]; col5 <- MAT[,5]
#col6 <- MAT[,6]; col7 <- MAT[,7]; col8 <- MAT[,8]; col9 <- MAT[,9]
#data <- data.frame(col1,col2,col3,col4,col5,col6,col7,col8,col9)
#data <- data.frame(MAT); data["X10"] = c(0,0,0,0,0,0,0,0)
data <- data.frame(MAT); data["X10"] = c(0,0,0,0)
names(data) <- c('NW Africa','NE Africa','Sahel','Middle East \n & C Asia','E Asia','N America','Australia','S America','S Africa','NA')

# regional emission budgets
par(mai=c(0.3,0.6,0.2,0.1), cex.axis=0.7)
Set2 = brewer.pal(8, 'Set2')[1:4]
# barplot with colors. Make sure that the plot and legends have same colors for items.
barplot(height=as.matrix(data)[,1:5], ylab="regional total emissions (Tg/yr)", beside=TRUE, col=Set2, ylim=c(0,800*3.7))
abline(h=0)
#legend(x=-0.5,y=880, legend = c("0.5x0.625-C19","0.9x1.25-C19","1.9x2.5-C19","4x5-C19","0.5x0.625-K14","0.9x1.25-K14","1.9x2.5-K14","4x5-K14"), cex=0.9, bty="n", fill=Set2, y.intersp=0.7)  #Add legends
legend(x=-0.5,y=(800*3.7), legend = c("0.5x0.625","0.9x1.25","1.9x2.5","4x5"), cex=0.9, bty="n", fill=Set2, y.intersp=0.9)  #Add legends
barplot(height=as.matrix(data)[,6:10], yaxt="n", beside=TRUE, col=Set2, ylim=c(0,800*3.7))
abline(h=0)

# global emission budgets
#MAT = cbind(sum(F_d.M2.int*area.M2/1e9),sum(F_d.C1.int*area.C1/1e9),sum(F_d.CESM.int*area.CESM/1e9),sum(F_d.GC.int*area.GC/1e9),sum(F_d.M2.noint*area.M2/1e9),sum(F_d.C1.noint*area.C1/1e9),sum(F_d.CESM.noint*area.CESM/1e9),sum(F_d.GC.noint*area.GC/1e9))
MAT = cbind(sum(F_d.M2.int*area.M2/1e9),sum(F_d.C1.int*area.C1/1e9),sum(F_d.CESM.int*area.CESM/1e9),sum(F_d.GC.int*area.GC/1e9)) * 4.5
data <- data.frame(MAT)
names(data) <- c("0.5x0.625","0.9x1.25","1.9x2.5","4x5")
#names(data) <- c('NW Africa','NE Africa','Sahel','Middle East','E Asia','N America','Australia','S America','S Africa')
#names(data) <- c("0.5x0.625-C19","0.9x1.25-C19","1.9x2.5-C19","4x5-C19","0.5x0.625-K14","0.9x1.25-K14","1.9x2.5-K14","4x5-K14")

par(mai=c(0.3,0.6,0.2,0.1))
Set2 = c(brewer.pal(8, 'Set2'))[1:4]
# barplot with colors. Make sure that the plot and legends have same colors for items.
barplot(height=as.matrix(data), ylab="global total emissions (Tg/yr)", beside=TRUE, col=Set2, ylim=c(0,2600*4.5), space=c(0.1,0.1,0.1,0.1), xaxt='n')
abline(h=0)
#legend(x=6,y=(2500*3.38), legend = c("0.5x0.625-C19","0.9x1.25-C19","1.9x2.5-C19","4x5-C19","0.5x0.625-K14","0.9x1.25-K14","1.9x2.5-K14","4x5-K14"), cex=1, bty="n", fill=Set2, y.intersp=0.8)  #Add legends
legend(x=3,y=(2600*4.5), legend = c("0.5x0.625","0.9x1.25","1.9x2.5","4x5"), cex=1, bty="n", fill=Set2, y.intersp=0.8)  #Add legends

#barplot(height=as.matrix(data)[,6:10], yaxt="n", beside=TRUE, col=Set2, ylim=c(0,800))
#abline(h=0)


###################################################################
# grid-by-grid emissions ratios between two different emissions

# 0.5 deg emission : 2 deg emission ratio
plot.field.log(F_d.M2.int[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)#; frame.region(frame='source', lwd=1.5)
F_d.M2.int.f19 = sp.dissolve(F_d.M2.int, lon.M2, lat.M2, LON.CESM, LAT.CESM)
plot.field.log(F_d.M2.int.f19[indlon,], LON.r, LAT.CESM, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)#; frame.region(frame='source', lwd=1.5)
plot.field.log(F_d.CESM.int[indlon,], LON.r, LAT.CESM, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)#; frame.region(frame='source', lwd=1.5)
A = F_d.M2.int.f19/F_d.CESM.int; A[which(A>1000)] = NaN
fct = sum(F_d.M2.int.f19,na.rm=T) / sum(F_d.CESM.int,na.rm=T)
plot.field.log((A/fct)[indlon,], LON.r, LAT.CESM, type='def', zlim=c(-1,2), legend.mar=5, col=TEMP_DIFF_65[c(seq(8,28,by=2),33,38:59)])#; frame.region(frame='source', lwd=1.5)

# 0.5 deg emission : 4 deg emission ratio
F_d.M2.int.f4 = sp.dissolve(F_d.M2.int, lon.M2, lat.M2, LON.GC, LAT.GC)
plot.field.log(F_d.M2.int.f4[i.GC,], LON.GC.r, LAT.GC, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)#; frame.region(frame='source', lwd=1.5)
plot.field.log(F_d.GC.int[i.GC,], LON.GC.r, LAT.GC, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)#; frame.region(frame='source', lwd=1.5)
A = F_d.M2.int.f4/F_d.GC.int; A[which(A>1000)] = NaN
fct = sum(F_d.M2.int.f4,na.rm=T) / sum(F_d.GC.int,na.rm=T)
plot.field.log((A/fct)[i.GC,], LON.GC.r, LAT.GC, type='def', zlim=c(-1,2), legend.mar=5, col=TEMP_DIFF_65[c(seq(8,28,by=2),33,38:59)])#; frame.region(frame='source', lwd=1.5)


# 0.5 deg emission : 1 deg emission ratio
#plot.field.log(F_d.M2.int, lon.M2, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)#; frame.region(frame='source', lwd=1.5)
F_d.M2.int.f09 = sp.dissolve(F_d.M2.int, lon.M2, lat.M2, LON.C1, LAT.C1)
plot.field.log(F_d.M2.int.f09[i.C1,], LON.C1.r, LAT.C1, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)#; frame.region(frame='source', lwd=1.5)
plot.field.log(F_d.C1.int[i.C1,], LON.C1.r, LAT.C1, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)#; frame.region(frame='source', lwd=1.5)
A = F_d.M2.int.f09/F_d.C1.int; A[which(A>1000)] = NaN
#A[which(A>5)] = 5
A[which(A>100)] = 100
fct = sum(F_d.M2.int.f09,na.rm=T) / sum(F_d.C1.int,na.rm=T)
#plot.field.log(A/fct, LON.C1, LAT.C1, type='def', zlim=c(-1,1), legend.mar=5, col=TEMP_DIFF_65)#; frame.region(frame='source', lwd=1.5)
plot.field.log((A/fct)[i.C1,], LON.C1.r, LAT.C1, type='def', zlim=c(-1,2), legend.mar=5, col=TEMP_DIFF_65[c(seq(8,28,by=2),33,38:59)])
plot.field.log((A/fct)[i.C1,], LON.C1.r, LAT.C1, type='def', zlim=c(-1,2), legend.mar=5, col=TEMP_DIFF_65[c(seq(8,28,by=2),33,38:59)], legend.only=T)
plot.field.proj(log10(A/fct)[i.C1,], LON.C1.r, LAT.C1, def.zlim=T, zlim=c(-1,2), col=TEMP_DIFF_65[c(seq(8,28,by=2),33,38:59)])

#####################################################
# save as .ncdf4 file

plot.field.log(A/fct, LON.C1, LAT.C1,type='def',zlim=c(-1,1),legend.mar=5) # ratio for rescaling emissions
plot.field.log(A/fct, LON.C1, LAT.C1,type='def',zlim=c(-1,1),legend.mar=5,col=TEMP_DIFF_19) # ratio
#plot.field.log(GWC.Ratio,LON,LAT,type='def',zlim=c(-5,5),legend.mar=5) # ratio

A = A/fct
A[which(is.na(A))]=1
plot.field.log(A, LON.C1, LAT.C1,type='def',zlim=c(-1,1),legend.mar=5) # ratio for rescaling emissions

files.clm = Sys.glob('/Volumes/GoogleDrive/My Drive/CESM_tempfiles/F2000climo/nemd5_F_L32.clm2.h1.2006-*.nc'); nc = nc_open(files.clm[1])
lon = ncvar_get(nc,'lon'); lat = ncvar_get(nc,'lat') # define dimensions

# array to be saved
arr.save = array(1, dim=c(length(lon), length(lat)))
arr.save[c(145:288,1:144),((2*18):(2*93))] = A            # global adjustment
plot.field.log(A,LON.C1,LAT.C1,type='def',zlim=c(-1,1),legend.mar=5,col=TEMP_DIFF_19) # ratio
plot.field.log(arr.save,lon,lat,type='def',zlim=c(-1,1),legend.mar=5,col=TEMP_DIFF_19,Pacific.centric=T) # ratio

# make sure the array to be saved is in a reasonable value such as from 0.01 to 10
quantile(arr.save,na.rm=T)

londim = ncdim_def('lon', 'degrees_east', as.numeric(lon), longname = 'Longitude')
latdim = ncdim_def('lat', 'degrees_north', as.numeric(lat), longname = 'Latitude')
#timedim = ncdim_def('date', 'YYYYMMDD', vals = as.numeric(date_vec), unlim = T, longname = 'Date vector')

# define variables
fillvalue = 1e32
dlname = c('emission rescaling factor due to resolution-dependence')
dummy1 = ncvar_def('DSTFLXT_CORRECTION', 'dimensionless', list(londim, latdim), fillvalue, dlname[1], prec='double')
#newfile = nc_create(filename='GWC_correction-CN-09x125-12222021.nc', var=list(dummy1), force_v4=T, verbose=T)  # when using this one, comment the others
newfile = nc_create(filename='DSTFLXT_correction-Globe_lim20-09x125-01192022.nc', var=list(dummy1), force_v4=T, verbose=T)  # when using this one, comment the others
ncvar_put(newfile, dummy1, arr.save) # put variables

# put additional attributes into dimensions and data variables
ncatt_put(newfile, 'lon', 'axis', 'X', prec='double')
ncatt_put(newfile, 'lat', 'axis', 'Y', prec='double')
#ncatt_put(newfile, 'date', 'axis', 't', prec='float')

# close the file, write data to disk
nc_close(newfile)



#nc = nc_open('GWC_correction-CN-09x125-12222021.nc')
nc = nc_open('DSTFLXT_correction-Globe-09x125-01162022.nc')
#nc = nc_open('DSTFLXT_correction-Globe-09x125-01192022.nc')
lonn = ncvar_get(nc,'lon')
latt = ncvar_get(nc,'lat')
DSTFLXT_CORRECTION = ncvar_get(nc,'DSTFLXT_CORRECTION')
plot.field(DSTFLXT_CORRECTION,lonn,latt,Pacific.centric=T,type='def',zlim=c(0,30))
plot.field.log(DSTFLXT_CORRECTION,lonn,latt,type='def',zlim=c(-1,1),legend.mar=5,col=TEMP_DIFF_19,Pacific.centric=T) # ratio

arr.save = DSTFLXT_CORRECTION
arr.save[which(arr.save>20)] = 20


# M2 vs CESM emission
plot.field.log(F_d.M2.int[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)#; frame.region(frame='source', lwd=1.5)
F_d.M2.int.f19 = sp.dissolve(F_d.M2.int, lon.M2, lat.M2, LON.CESM, LAT.CESM)
plot.field.log(F_d.CESM.int, LON.CESM.r, LAT.CESM, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)  
plot.field.log(F_d.M2.int.f19/F_d.CESM.int, LON.CESM, LAT.CESM, type='def', zlim=c(-1,1), legend.mar=5, col=WBGYR)

par(mai=c(0.6,0.6,0.2,0.2))
plot(NA,NA, xlab='1.9°x2.5° emissions (kg/m2/yr)', ylab='0.5°x0.625° emissions (kg/m2/yr)',xlim=c(1e-6,1), ylim=c(1e-6,1), xaxt='n', yaxt='n',log='xy')
axis(side = 1, at=(10^((-6):0)), labels=expression('10'^-6,'10'^-5,'10'^-4,'10'^-3,'10'^-2,'10'^-1,'0')); axis(side = 2, at=(10^((-6):0)), labels=expression('10'^-6,'10'^-5,'10'^-4,'10'^-3,'10'^-2,'10'^-1,'0'))
abline(a=0,b=1, lwd=2); abline(a=-1,b=1, lwd=2, lty=2); abline(a=-2,b=1, lwd=2, lty=3); abline(a=1,b=1, lwd=2, lty=2); abline(a=2, b=1, lwd=2, lty=3)
#abline(v=-4, col='grey', lwd=2); abline(h=-4, col='grey', lwd=2)
points((F_d.CESM.int), (F_d.M2.int.f19), 'p', col=TEMP_DIFF_65[6], pch=19)
#points((F_d.CESM.int), (F_d.M2.int.f19), 'p', col='black',log='xy', xlim=c(1e-6,1), ylim=c(1e-6,1))


plot(log10(F_d.CESM.int), log10(F_d.M2.int.f19), 'p', xlim=c(-6,0), ylim=c(-6,0), col='blue', xlab='1.9x2.5 emission (kg/m2/yr)', ylab='0.5x0.625 emission (kg/m2/yr)', log='x')


plot.field.log(F_d.M2.int.f19/F_d.CESM.int, LON.CESM, LAT.CESM, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR, horizontal = T)


#########################################################

# correction map seasonality


d.MAM = make.date.vec(20060301, 20060531)
FF = get.tot.dustemis(d.MAM, dt, res='05x0625', use.intermittency='TRUE'); 
F_d.M2.MAM = FF$F_d
sum(F_d.M2.MAM*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
FF = get.tot.dustemis(d.MAM, dt, res='1x1', use.intermittency=TRUE); 
F_d.C1.MAM = FF$F_d
sum(F_d.C1.MAM*area.C1, na.rm=TRUE)/1e9   # Tg/yr-1

d.JJA = make.date.vec(20060601, 20060831)
FF = get.tot.dustemis(d.JJA, dt, res='05x0625', use.intermittency='TRUE'); 
F_d.M2.JJA = FF$F_d
sum(F_d.M2.JJA*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
FF = get.tot.dustemis(d.JJA, dt, res='1x1', use.intermittency=TRUE); 
F_d.C1.JJA = FF$F_d
sum(F_d.C1.JJA*area.C1, na.rm=TRUE)/1e9   # Tg/yr-1

d.SON = make.date.vec(20060901, 20061130)
FF = get.tot.dustemis(d.SON, dt, res='05x0625', use.intermittency='TRUE'); 
F_d.M2.SON = FF$F_d
sum(F_d.M2.SON*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
FF = get.tot.dustemis(d.SON, dt, res='1x1', use.intermittency=TRUE); 
F_d.C1.SON = FF$F_d
sum(F_d.C1.SON*area.C1, na.rm=TRUE)/1e9   # Tg/yr-1

d.DJF = c(make.date.vec(20060101, 20060228),make.date.vec(20061201, 20061231))
FF = get.tot.dustemis(d.DJF, dt, res='05x0625', use.intermittency='TRUE'); 
F_d.M2.DJF = FF$F_d
sum(F_d.M2.DJF*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
FF = get.tot.dustemis(d.DJF, dt, res='1x1', use.intermittency=TRUE); 
F_d.C1.DJF = FF$F_d
sum(F_d.C1.DJF*area.C1, na.rm=TRUE)/1e9   # Tg/yr-1

# MAM correction 
F_d.M2.MAM.f09 = sp.dissolve(F_d.M2.MAM, lon.M2, lat.M2, LON.C1, LAT.C1)
A = F_d.M2.MAM.f09/F_d.C1.MAM; A[which(A>1000)] = NaN; A[which(F_d.M2.MAM.f09<0.25e-4)] = NaN; #
fct = sum(F_d.M2.MAM.f09,na.rm=T) / sum(F_d.C1.MAM,na.rm=T)
plot.field.log((A/fct)[i.C1,], LON.C1.r, LAT.C1, type='def', zlim=c(-1,2), legend.mar=5, col=TEMP_DIFF_65[c(seq(8,28,by=2),33,38:59)])#; frame.region(frame

# JJA correction
F_d.M2.JJA.f09 = sp.dissolve(F_d.M2.JJA, lon.M2, lat.M2, LON.C1, LAT.C1)
A = F_d.M2.JJA.f09/F_d.C1.JJA; A[which(A>1000)] = NaN; #A[which(F_d.M2.JJA.f09<0.25e-4)] = NaN; #
fct = sum(F_d.M2.JJA.f09,na.rm=T) / sum(F_d.C1.JJA,na.rm=T)
plot.field.log((A/fct)[i.C1,], LON.C1.r, LAT.C1, type='def', zlim=c(-1,2), legend.mar=5, col=TEMP_DIFF_65[c(seq(8,28,by=2),33,38:59)])#; frame.region(frame

# SON correction
F_d.M2.SON.f09 = sp.dissolve(F_d.M2.SON, lon.M2, lat.M2, LON.C1, LAT.C1)
A = F_d.M2.SON.f09/F_d.C1.SON; A[which(A>1000)] = NaN; #A[which(F_d.M2.SON.f09<0.25e-4)] = NaN; 
fct = sum(F_d.M2.SON.f09,na.rm=T) / sum(F_d.C1.SON,na.rm=T)
plot.field.log((A/fct)[i.C1,], LON.C1.r, LAT.C1, type='def', zlim=c(-1,2), legend.mar=5, col=TEMP_DIFF_65[c(seq(8,28,by=2),33,38:59)])#; frame.region(frame

# DJF correction
F_d.M2.DJF.f09 = sp.dissolve(F_d.M2.DJF, lon.M2, lat.M2, LON.C1, LAT.C1)
A = F_d.M2.DJF.f09/F_d.C1.DJF; A[which(A>1000)] = NaN; #A[which(F_d.M2.DJF.f09<0.25e-4)] = NaN #
fct = sum(F_d.M2.DJF.f09,na.rm=T) / sum(F_d.C1.DJF,na.rm=T)
plot.field.log((A/fct)[i.C1,], LON.C1.r, LAT.C1, type='def', zlim=c(-1,2), legend.mar=5, col=TEMP_DIFF_65[c(seq(8,28,by=2),33,38:59)])#; frame.region(frame

# Annual correction
A = (F_d.M2.DJF.f09+F_d.M2.MAM.f09+F_d.M2.JJA.f09+F_d.M2.SON.f09)/(F_d.C1.DJF+F_d.C1.MAM+F_d.C1.JJA+F_d.C1.SON); A[which(A>1000)] = NaN#; A[which((F_d.M2.DJF.f09+F_d.M2.MAM.f09+F_d.M2.JJA.f09+F_d.M2.SON.f09)<1e-4)] = NaN #
fct = sum((F_d.M2.DJF.f09+F_d.M2.MAM.f09+F_d.M2.JJA.f09+F_d.M2.SON.f09),na.rm=T) / sum((F_d.C1.DJF+F_d.C1.MAM+F_d.C1.JJA+F_d.C1.SON),na.rm=T)
plot.field.log((A/fct)[i.C1,], LON.C1.r, LAT.C1, type='def', zlim=c(-1,2), legend.mar=5, col=TEMP_DIFF_65[c(seq(8,28,by=2),33,38:59)])#; frame.region(frame

plot.field.log((A/fct)[i.C1,], LON.C1.r, LAT.C1, type='def', zlim=c(-1,2), legend.mar=5, col=TEMP_DIFF_65[c(seq(8,28,by=2),33,38:59)], horizontal=T)#; frame.region(frame
