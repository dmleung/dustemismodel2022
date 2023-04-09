#############################################################
#############################################################

# 5 Oct 2021
# dust emission code analysis
# functions for plot and analysis

# 18 Nov 2021
# The code has been revised to include Zender et al. (2003) scheme as the dust production model (DPM). Other codes mostly only have Kok et al. (2014) scheme as the DPM.

# set working directory
setwd('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR')
# load libraries 
library(ncdf4); library(fields); library(maps); library(Metrics); library(pracma); library(abind)
# load functions
source("get_geo.R"); source("get_met.R"); source("sptial_plot_fns.R")
source("dust_research_fns.R")
# load color scheme
load("WBGYR_scheme.RData"); load('TEMP_DIFF.RData'); load('LEAF_scheme.RData'); load('CBR_DRYWET.RData'); load('BAGYOR_scheme.RData')
# load landfilt 
.env = new.env()
load('landfilt_05x0625.RData',envir=.env)
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
# redefine landfilt coordinates
landfilt.M2 = .env$landfilt.M2[indi,indj]
# choose land cover averaging method: LC0 (rock and veg) or LC1 (rock, veg and mix)
LC = 'LC0'
##################
# calculate model parameters based on user inputs
date_vec = make.date.vec(startdate, enddate)  # make date vec
time = 1:(length(date_vec)*24/dt)  # a vector of timesteps
indt = seq(dt, 24, by=dt)   # hour of day
# get lon lat
nc = nc_open('/Volumes/GoogleDrive/My Drive//MERRA2/upscale_2006/MERRA2_300.tavg1_2d_lnd_Nx.20060101.SUB.nc')
lon.M2 = ncvar_get(nc, 'lon')[indi]; lat.M2 = ncvar_get(nc, 'lat')[indj]
nc_close(nc)

i.M2 = c(20:576,1:19)
lon.M2.r = lon.M2[i.M2]; lon.M2.r[558:576]=lon.M2.r[558:576]+360

##################
# get a 2-D area map for MERRA2 resolution
dlon = lon.M2[2]-lon.M2[1]; dlat = lat.M2[2]-lat.M2[1]
lat.area = NULL
for (j in 1:length(lat.M2)) {
  lat1 = lat.M2[j]-dlat/2; lat2 = lat.M2[j]+dlat/2
  lon1 = lon.M2[j]-dlon/2; lon2 = lon.M2[j]+dlon/2
  lat.area = c(lat.area, area.latlon(lat1, lon1, lat2, lon2))
}
area.M2 = NULL; for (i in 1:length(lon.M2)) {area.M2 = rbind(area.M2, lat.area)}
area.M2 = area.M2*1e6
#plot.field(area.M2, lon.M2, lat.M2, legend.mar=5)


##################
# function for getting annual total dust emissions
get.tot.dustemis = function(date_vec, dt, directory) {
  # get annual dust emission total per grid (kg m-2 yr-1) 
  .env1 = new.env()
  for (d in 1:length(date_vec)) {
    print(date_vec[d])
    #filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')
    #filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006_new/', directory, '/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')
    filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006_new/', directory, '/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')
    load(filename, envir=.env1)
    #if (d==1) {lon.M2 = .env1$lon.M2; lat.M2 = .env1$lat.M2; F_d = .env1$F_d}
    #else if (d >= 2) {F_d = F_d + .env1$F_d}
    if (d==1) {F_d=.env1$F_d; F_d[which(is.na(F_d))]=0}
    else if (d >= 2) {F_d1=.env1$F_d; F_d1[which(is.na(F_d1))]=0; F_d=F_d+F_d1}
    #F_d = abind(F_d, (apply(.env1$F_d, c(1,2), sum, na.rm=T)*3600*dt), along=3)
    #F_d = apply(F_d, c(1,2), sum, na.rm=T)
  }
  dim(F_d)
  F_d = apply(F_d, c(1,2), sum, na.rm=T)*3600*dt
  return(list(F_d=F_d,lon.M2=lon.M2, lat.M2=lat.M2))
  print(paste('total global dust emission is ',sum(F_d*area.M2, na.rm=TRUE)/1e9, ' Tg/yr',sep=''))
}

FF = get.tot.dustemis(date_vec, dt=2, directory='K14_default')
F_d.K14_def = FF$F_d
FF = get.tot.dustemis(date_vec, dt=2, directory='K14_u_ft_Dp=130')
F_d.K14_130 = FF$F_d
FF = get.tot.dustemis(date_vec, dt=2, directory='K14_u_ft_Dp=130_F_eff')
F_d.K14_uft_F_eff = FF$F_d
FF = get.tot.dustemis(date_vec, dt=2, directory='K14_u_it_Dp=130_F_eff')
F_d.K14_F_eff = FF$F_d
FF = get.tot.dustemis(date_vec, dt=2, directory='L21_final')
F_d.L21 = FF$F_d
save('F_d.K14_def','F_d.K14_130','F_d.K14_uft_F_eff','F_d.K14_F_eff','F_d.L21','lon.M2.r','lat.M2','date_vec',file='MERRA2_dustmod_emissions_05x06_c16Jan22.RData')

# a) K14 emissions
plot.field.log(F_d.K14_def[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=5, col=WBGYR, ps=15)
sum(F_d.K14_def*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
#plot.field.log(F_d.K14_def[i.M2,]/29254*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

# b) L21 emissions
plot.field.log(F_d.L21[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=5, col=WBGYR, ps=15)
sum(F_d.L21*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
#plot.field.log(F_d.L21[i.M2,]/11748*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

# c) absolute change (unnormalized) from a) from Dp = 127um
plot.field((F_d.K14_130 - F_d.K14_def)[i.M2,], lon.M2.r, lat.M2, type='sign', legend.mar=5, col=TEMP_DIFF_65, zlim=c(-2,2), ps=16)
sum(F_d.K14_130*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log.diff((F_d.K14_130)[i.M2,], (F_d.K14_def)[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-4)
#plot.field.log(-(F_d.K14_130 - F_d.K14_def)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=4.5, col=TEMP_DIFF_65[seq(33,5,by=(-1))])

# d) relative change from a) (ratio) after normalization from Dp = 127um 
R1 = (F_d.K14_130/23984 / (F_d.K14_def/29254))
R1[which(R1 > 1000)] = NaN
#R1[which(R1<0.8)] = R1[which(R1<0.8)]*1.4
.col = TEMP_DIFF_65[c(9:31,33,35:57)]
plot.field.log(R1[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-1,1), legend.mar=5, col=.col, ps=15)
plot.field.log((F_d.K14_130/F_d.K14_def)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-1,1), legend.mar=5, col=TEMP_DIFF_65, ps=15)

# e) absolute change (unnormalized) from c) from including Feff
plot.field((F_d.K14_uft_F_eff - F_d.K14_130)[i.M2,], lon.M2.r, lat.M2, type='sign', legend.mar=5, col=TEMP_DIFF_65, zlim=c(-2,2), ps=15)
plot.field.log.diff((F_d.K14_uft_F_eff)[i.M2,], (F_d.K14_130)[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-4, oom.max=10)
#plot.field.log(-(F_d.K14_uft_F_eff - F_d.K14_130)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=4.5, col=TEMP_DIFF_65[seq(33,5,by=(-1))])

# f) relative change (ratio) from c) after normalization from including Feff
R2 = (F_d.K14_uft_F_eff/2886 / (F_d.K14_130/23984))
sum(F_d.K14_uft_F_eff*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
#R2[which(F_d.K14_uft_F_eff/2884*5000 < 0.0001)] = NaN
#R2[which(R2 > 1000)] = NaN
.col = TEMP_DIFF_65[c(7:33,seq(35,59,by=2))]
plot.field.log(R2[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-2,1), legend.mar=5, col=.col, ps=15)
plot.field.log((F_d.K14_uft_F_eff/F_d.K14_130)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-1,1), legend.mar=5, col=TEMP_DIFF_65, ps=15)


# absolute change (unnormalized) from using impact threshold
plot.field.log.diff((F_d.K14_F_eff)[i.M2,], (F_d.K14_uft_F_eff)[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-4, oom.max=10)

# absolute change (unnormalized) from including intermittency
plot.field.log.diff((F_d.L21)[i.M2,], (F_d.K14_F_eff)[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-4, oom.max=0.1)

# g) absolute change (unnormalized) from e) from including intermittency
plot.field((F_d.L21 - F_d.K14_uft_F_eff)[i.M2,], lon.M2.r, lat.M2, type='sign', legend.mar=5, col=TEMP_DIFF_65, zlim=c(-2,2), ps=15)
plot.field.log.diff((F_d.L21)[i.M2,], (F_d.K14_uft_F_eff)[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-4, oom.max=10)
#plot.field.log((F_d.L21 - F_d.K14_uft_F_eff)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=5, col=TEMP_DIFF_65[33:61])

# h) relative change (ratio) from e) after normalization from including intermittency
sum(F_d.K14_F_eff*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
sum(F_d.L21*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
R5 = (F_d.L21/11748 / (F_d.K14_uft_F_eff/2886))
R5[which(F_d.L21/11748*5000 < 1e-4)] = NaN
#R3[which(R3 > 1e50)] = NaN
.col = TEMP_DIFF_65[c(seq(11,33,by=2), 35:56)]
plot.field.log(R5[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-1,2), legend.mar=5, col=.col)
.col = TEMP_DIFF_65[c(seq(11,33,by=3),33, 35:57)]
plot.field.log(R5[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-1,3), legend.mar=5, col=.col)
.col = TEMP_DIFF_65[c(seq(11,31,by=4),33,35:47 ,49,51:60)]
plot.field.log(R5[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-1,4), legend.mar=5, col=.col, ps=14)

#i) absolute change from a) to b) 
plot.field((F_d.L21 - F_d.K14_def)[i.M2,], lon.M2.r, lat.M2, type='sign', legend.mar=5, col=TEMP_DIFF_65, zlim=c(-2,2), ps=15)

plot.field.log.diff((F_d.L21)[i.M2,], (F_d.K14_def)[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-4, oom.max=10)

A = (F_d.L21 - F_d.K14_def)
#A[which(A<0)] = -1; A[which(A>=0)] = 1
A[which(A<0)] = 0
#plot.field.log(A[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=4.5, col=TEMP_DIFF_65[33:61])
A = log(A); A[which(is.infinite(A))] = NaN
plot.field(A[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-5,1), legend.mar=4.5, col=TEMP_DIFF_65[33:61])
range(A,na.rm=T)
B = (F_d.K14_def - F_d.L21); B[which(B<0)] = 0
B = log(B); B[which(is.infinite(B))] = NaN
range(B,na.rm=T)
plot.field(B[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-5,1), legend.mar=4.5, col=TEMP_DIFF_65[seq(33,5,by=(-1))])

plot.field(B[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-5,1), legend.mar=4.5, col=TEMP_DIFF_65[seq(33,5,by=(-1))])
#plot.field((log(A-B))[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=4.5, col=WBGYR)
C = A+4.691361; C[which(C<0)] = 0
B1 = B+4.691361; B1[which(B1<0)] = 0
C[which(!is.na(B1))] = (-B1)[which(!is.na(B1))]
range(C,na.rm=T)
plot.field(C[i.M2,], lon.M2.r, lat.M2, legend.mar=4.5, col=TEMP_DIFF_65[5:61], type='sign')

#j) relative change (ratio) after normalization from a) to b)
R6 = (F_d.L21/11748 / (F_d.K14_def/29254))
R6[which(F_d.L21/11748*5000 < 1e-4)] = NaN
#R3[which(R3 > 1e50)] = NaN
.col = TEMP_DIFF_65[c(seq(9,31,by=2),33, 35:58)]
plot.field.log(R6[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-1,2), legend.mar=5, col=.col, ps=15)

#.col = TEMP_DIFF_65[c(seq(11,31,by=4),33,35:47,49,51:60)]
.col = TEMP_DIFF_65[c(seq(7,31,by=2),33,35,37,39:47,49,51:58)]
plot.field.log(R6[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-2,3), legend.mar=5, col=.col, ps=14)

# Figure S6
# a) K14
plot.field.log(F_d.K14_def[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-4,1), legend.mar=5, col=WBGYR, ps=15)
sum(F_d.K14_def*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1

# b) 
plot.field.log(F_d.K14_130[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-4,1), legend.mar=5, col=WBGYR, ps=15)
sum(F_d.K14_130*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1

# c) 
plot.field.log(F_d.K14_uft_F_eff[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-4,1), legend.mar=5, col=WBGYR, ps=15)
sum(F_d.K14_uft_F_eff*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1

# d) 
plot.field.log(F_d.K14_F_eff[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-4,1), legend.mar=5, col=WBGYR, ps=15)
sum(F_d.K14_F_eff*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1

# e) 
plot.field.log(F_d.L21[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-4,1), legend.mar=5, col=WBGYR, ps=15)
sum(F_d.L21*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1

#FF = get.tot.dustemis(date_vec, dt=2, directory='K14_default')
#F_d.K14_def = FF$F_d
sum(F_d.K14_def*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.K14_def[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=5, col=WBGYR)
plot.field.log(F_d.K14_def[i.M2,]/29254*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
plot.field.log(F_d.K14_def[i.M2,]/29254*5000, lon.M2.r, lat.M2, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR)
#sum(F_d*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1

FF = get.tot.dustemis(date_vec, dt=2, directory='K14_u_ft_Dp=130')
F_d.K14_130 = FF$F_d
sum(F_d.K14_130*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
#plot.field.log(F_d.K14_130[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=5)
plot.field.log(F_d.K14_130[i.M2,]/23984*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
R1 = (F_d.K14_130/23984 / (F_d.K14_def/29254))
R1[which(R1 > 1000)] = NaN
#R1[which(R1<0.8)] = R1[which(R1<0.8)]*1.4
.col = TEMP_DIFF_65[c(9:30,33,36:57)]
plot.field.log(R1[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-1,1), legend.mar=5, col=.col)

FF = get.tot.dustemis(date_vec, dt=2, directory='K14_u_ft_Dp=130_F_eff')
F_d.K14_uft_F_eff = FF$F_d
sum(F_d.K14_uft_F_eff*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.K14_uft_F_eff[i.M2,]/2886*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
R2 = (F_d.K14_uft_F_eff/2886 / (F_d.K14_130/23984))
#R2[which(F_d.K14_uft_F_eff/2884*5000 < 0.0001)] = NaN
#R2[which(R2 > 1000)] = NaN
.col = TEMP_DIFF_65[c(7:33,seq(35,59,by=2))]
plot.field.log(R2[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-2,1), legend.mar=5, col=.col)

#FF = get.tot.dustemis(date_vec, dt=2, directory='K14_u_it_Dp=130')
#F_d.K14_uit = FF$F_d
#sum(F_d.K14_uit*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
#plot.field.log(F_d.K14_uit[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=5)
#plot.field.log(F_d.K14_uit[i.M2,]/87578*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

FF = get.tot.dustemis(date_vec, dt=2, directory='K14_u_it_Dp=130_F_eff')
F_d.K14_F_eff = FF$F_d
sum(F_d.K14_F_eff*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
#plot.field.log(F_d.K14_F_eff[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=5)
plot.field.log(F_d.K14_F_eff[i.M2,]/13181*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
R3 = (F_d.K14_F_eff/13181 / (F_d.K14_uft_F_eff/2886))
R3[which(F_d.K14_F_eff/13181*5000 < 0.0001)] = NaN
#R3[which(R3 > 1e50)] = NaN
.col = TEMP_DIFF_65[c(seq(9,33,by=2), 34:57)]
plot.field.log(R3[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-1,2), legend.mar=5, col=.col)

FF = get.tot.dustemis(date_vec, dt=2, directory='L21_final')
F_d.L21 = FF$F_d
#plot.field.log(F_d.L21[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-3,1), legend.mar=5, col=WBGYR)
sum(F_d.L21*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.L21[i.M2,]/11748*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
plot.field.log(F_d.L21[i.M2,]/11748*5000, lon.M2.r, lat.M2, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR)
plot.field.log(F_d.L21[i.M2,]/11748*5000, lon.M2.r, lat.M2, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR[c(1,16,30,45,60,70,90,150,180,220)])
R4 = (F_d.L21/11748 / (F_d.K14_F_eff/13181))
#R3[which(F_d.K14_F_eff/11828*5000 < 0.001)] = NaN
#R4[which((R4 > 1e20))] = NaN
.col = TEMP_DIFF_65[c(11:33,seq(35,55,by=2))]
plot.field.log(R4[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-2,1), legend.mar=5, col=.col)


# with and without C19 scheme
R5 = (F_d.L21/11748 / (F_d.K14_uft_F_eff/2886))
R5[which(F_d.L21/11748*5000 < 1e-4)] = NaN
#R3[which(R3 > 1e50)] = NaN
.col = TEMP_DIFF_65[c(seq(9,31,by=2), 36:57)]
plot.field.log(R5[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-1,2), legend.mar=5, col=.col)


# with and without all of our modifications
R6 = (F_d.L21/11748 / (F_d.K14_def/29254))
R6[which(F_d.L21/11748*5000 < 1e-4)] = NaN
#R3[which(R3 > 1e50)] = NaN
.col = TEMP_DIFF_65[c(seq(9,31,by=2), 36:57)]
plot.field.log(R6[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-1,2), legend.mar=5, col=.col)




F_d.int.reg = extract.value.region(F_d.int*area.M2/1e9, lon=lon.M2, lat=lat.M2, metric='sum')
F_d.noint.reg = extract.value.region(F_d.noint*area.M2/1e9, lon=lon.M2, lat=lat.M2, metric='sum')

# plot emissions for Z03-LC0 scheme 
FF = get.tot.dustemis(date_vec, dt, DPM='Z03')
#lon.M2 = FF$lon.M2; lat.M2 = FF$lat.M2; 
F_d.Z03.G = FF$F_d
plot.field.log(F_d.Z03.G, lon.M2, lat.M2, type='def', zlim=c(-4,0), legend.mar=5)
plot.field.log(F_d.Z03.G, lon.M2, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)


##################
# function for getting annual mean F_eff
get.F_eff = function(date_vec, dt) {
  # get annual dust emission total per grid (kg m-2 yr-1) 
  .env1 = new.env()
  for (d in 1:length(date_vec)) {
    print(date_vec[d])
    #filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_F_eff_',LC,'/F_eff_hybrid_05x0625_', date_vec[d], '.RData', sep='')
    filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006/OUTPUTS_F_eff_',LC,'/F_eff_hybrid_05x0625_', date_vec[d], '.RData', sep='')
    load(filename, envir=.env1)
    if (d==1) {lon.M2 = .env1$lon.M2; lat.M2 = .env1$lat.M2; F_eff = .env1$F_eff.LC}
    else if (d >= 2) {F_eff = F_eff + .env1$F_eff.LC}
    #F_d = abind(F_d, (apply(.env1$F_d, c(1,2), sum, na.rm=T)*3600*dt), along=3)
    #F_d = apply(F_d, c(1,2), sum, na.rm=T)
  }
  F_eff = F_eff/length(date_vec)
  return(list(F_eff=F_eff,lon.M2=lon.M2, lat.M2=lat.M2))
}

FF = get.F_eff(date_vec, dt)
lon.M2 = FF$lon.M2; lat.M2 = FF$lat.M2; F_eff = FF$F_eff
plot.field(F_eff[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1), col=WBGYR)
plot.field.log(F_d.noint, lon.M2, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

##################
# also can get and plot MERRA2 fields within specified time period
files_flx = files_slv = files_lnd = NULL   # character strings storing the file names
for (dd in 1:length(date_vec)) {
  date_vec[dd]
  files_flx = c(files_flx, Sys.glob(paste('/Volumes/GoogleDrive/My Drive/upscale_2006/tavg1_2d_flx_Nx/MERRA2*flx*', date_vec[dd], '*.nc', sep='')))  # for surface flux variables: surface USTAR, Z0M, and RHOA (air density)
  files_slv = c(files_slv, Sys.glob(paste('/Volumes/GoogleDrive/My Drive/upscale_2006/tavg1_2d_slv_Nx/MERRA2*slv*', date_vec[dd], '*.nc', sep='')))  # for surface level variables: surface U10, V10, and DISPH (displacement height)
  files_lnd = c(files_lnd, Sys.glob(paste('/Volumes/GoogleDrive/My Drive/upscale_2006/tavg1_2d_lnd_Nx/MERRA2*lnd*', date_vec[dd], '*.nc', sep='')))  # for land variables: SFMC (Volumetric water content) and LAI
}

# not sure why sometimes vectors are repeating the same file name twice
if (length(files_flx) != length(date_vec)) {
  which(duplicated(files_flx)); files_flx = files_flx[-which(duplicated(files_flx))]
  which(duplicated(files_slv)); files_slv = files_slv[-which(duplicated(files_slv))]
  which(duplicated(files_lnd)); files_lnd = files_lnd[-which(duplicated(files_lnd))]  
}

i1 = indi[1]; j1 = indj[1]; t1 = indt[1]
il = length(indi); jl = length(indj); tl = length(indt)

# function for getting annual mean MERRA2 fields
get.merra2.field = function(files, date_vec, dt, var=NULL) {
  # get annual dust emission total per grid (kg m-2 yr-1) 
  strings = strsplit(files, '[.]')
  date_files = as.double(lapply(strings, `[[`, 3)) # date
  
  ind = match(date_vec, date_files)
  print(ind)
  .env1 = new.env()
  for (d in 1:length(date_vec)) {
    print(date_vec[d])
    nc = nc_open(files[ind[d]])#; names(nc$var)#, start=c(11,26), count=c(5,5)
    
    #if (d==1) {VAR = ncvar_get(nc, var, start=c(i1,j1,t1), count=c(il,jl,tl))} # m/s
    #else if (d >= 2) {VAR = VAR + ncvar_get(nc, var, start=c(i1,j1,t1), count=c(il,jl,tl))}
    if (d==1) {VAR = ncvar_get(nc, var, start=c(i1,j1,1), count=c(il,jl,24))[,,indt]} # m/s
    else if (d >= 2) {VAR = VAR + ncvar_get(nc, var, start=c(i1,j1,1), count=c(il,jl,24))[,,indt]}
    
  }
  #dim(VAR)
  VAR = apply(VAR, c(1,2), mean, na.rm=T)/length(date_vec)
  return(list(VAR=VAR))
}

# get USTAR
date_vec1 = date_vec[-308]
FF = get.merra2.field(files_flx, date_vec1, dt, var='USTAR')
FF = get.merra2.field(files_flx, date_vec, dt, var='USTAR')
ustar = FF$VAR
plot.field((ustar*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,0.6))
plot.field((ustar*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,0.5), col=BAGYOR)
plot.field((ustar*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,0.5), col=WBGYR)

lon.M2.f = lon.M2[c(289:576,1:288)]; lon.M2.f[which(lon.M2.f<0)] = lon.M2.f[which(lon.M2.f<0)]+360
gwc_sfc.f09 = sp.dissolve(gwc_sfc[c(289:576,1:288),], lon.M2.f, lat.M2, LON.C1, LAT.C1)
landfrac.C1[which(!is.na(landfrac.C1))] = 1
plot.field.log(gwc_sfc.f09*landfrac.C1, LON.C1, LAT.C1, legend.mar=4, type='def', zlim=c(-2,0.6), col=tim.colors(32)[32:1], Pacific.centric = T)


# get LAI (seasonally)
d.MAM = make.date.vec(20060301, 20060531)
FFMAM = get.merra2.field(files_lnd, d.MAM, dt, var='LAI')
LAI.MAM = FFMAM$VAR
LAI.MAM1 = filt.art.data(LAI.MAM, q = 0.3, screen = 'low')  # take out the low values around the continents
plot.field((LAI.MAM1*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1.01), col=LEAF_256[c(3:210,220:222)])
plot.field((LAI.MAM1*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1.01), col=MPL_YG_128[c(3:115)])


d.JJA = make.date.vec(20060601, 20060831)
FFJJA = get.merra2.field(files_lnd, d.JJA, dt, var='LAI')
LAI.JJA = FFJJA$VAR
LAI.JJA1 = filt.art.data(LAI.JJA, q = 0.3, screen = 'low')
plot.field((LAI.JJA1*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1.01), col=LEAF_256[c(3:210,220:222)])
plot.field((LAI.JJA1*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1.01), col=MPL_YG_128[c(3:115)])

d.SON = make.date.vec(20060901, 20061130)
FFSON = get.merra2.field(files_lnd, d.SON, dt, var='LAI')
LAI.SON = FFSON$VAR
LAI.SON1 = filt.art.data(LAI.SON, q = 0.3, screen = 'low')
plot.field((LAI.SON1*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1.01), col=LEAF_256[c(3:210,220:222)])
plot.field((LAI.SON1*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1.01), col=MPL_YG_128[c(3:115)])

d.DJF = c(make.date.vec(20060101, 20060228), make.date.vec(20061201, 20061231))
FFDJF = get.merra2.field(files_lnd, d.DJF, dt, var='LAI')
LAI.DJF = FFDJF$VAR
LAI.DJF1 = filt.art.data(LAI.DJF, q = 0.3, screen = 'low')
plot.field((LAI.DJF1*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1.01), col=LEAF_256[c(3:210,220:222)])
plot.field((LAI.DJF1*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1.01), col=MPL_YG_128[c(3:115)])
#plot.field.log((LAI.DJF1*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(-1,0), col=MPL_YG_128[c(3:115)])

.col=two.colors(n=65, start=larry.colors()[1], end='darkgreen', middle = 'yellow',alpha = 0.5)

plot.field((LAI.DJF1*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1.01), col=LEAF_256[c(3:210,230:232)])
plot.field((LAI.DJF1*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1), col=MPL_YG_128[c(3:115)])



# annual mean LAI

FF = get.merra2.field(files_lnd, date_vec, dt, var='LAI')

LAI.ann = 0.25 * (LAI.DJF1 + LAI.MAM1 + LAI.JJA1 + LAI.SON1)
plot.field((LAI.ann*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1.01), col=LEAF_256[c(3:210,220:222)])
plot.field((LAI.ann*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1), col=LEAF_256[c(3:190)])
plot.field((LAI.ann*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,1.01), col=MPL_YG_128[c(3:115)])




# get soil moisture
date_vec = c(20060101,20060115,20060201,20060215,20060301,20060315,20060401,20060415,20060501,20060515,20060601,20060615,20060701,20060715,20060801,20060815,20060901,20060915,20061001,20061015,20061101,20061115,20061201,20061215)
FF = get.merra2.field(files_lnd, date_vec, dt, var='SFMC')
SFMC = FF$VAR
plot.field(SFMC[i.M2,]*1000*0.05/5*10, lon.M2.r, lat.M2, type='def', zlim=c(-2,48), legend.mar=3, col=grads_rainbow[1:12]) # m3water/m3soil to kgwater/m3soil * msoil
#plot.field((SFMC*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,0.5))
#plot.field.log((SFMC*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=5, type='def', zlim=c(-2.1,0))
plot.field((SFMC*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(0,0.5), col=tim.colors(32)[32:1])
plot.field.log((SFMC*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=4, type='def', zlim=c(-2,1))
plot.field.log((SFMC*landfilt.M2)[i.M2,], lon.M2.r, lat.M2, legend.mar=4, type='def', zlim=c(-2.1,0), col=CBR_DRYWET_65[c(1,seq(5,64,by=5),64)])
nc = nc_open('/Volumes/GoogleDrive/My Drive/MERRA2/upscale_2006/MERRA2_100.const_2d_lnd_Nx.00000000.nc4')#; names(nc$var)
poros = ncvar_get(nc, 'poros')[indi,indj] # porosity (m3 / m3, or volumetric soil water) 
plot.field(poros[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0.34,0.5))
dns_slt = 2650    # soil grain density
bulk_den = (1 - poros)*dns_slt 
SHR_CONST_RHOFW = 1000
gwc_sfc = SHR_CONST_RHOFW*sweep(SFMC, c(1,2), bulk_den, '/')
plot.field(gwc_sfc[i.M2,], lon.M2.r, lat.M2, legend.mar=3, type='def', zlim=c(0,0.35))
plot.field.log(gwc_sfc[i.M2,], lon.M2.r, lat.M2, legend.mar=4, type='def', zlim=c(-2,0.6), col=tim.colors(32)[32:1])
lon.M2.f = lon.M2[c(289:576,1:288)]; lon.M2.f[which(lon.M2.f<0)] = lon.M2.f[which(lon.M2.f<0)]+360
gwc_sfc.f09 = sp.dissolve(gwc_sfc[c(289:576,1:288),], lon.M2.f, lat.M2, LON.C1, LAT.C1)
landfrac.C1[which(!is.na(landfrac.C1))] = 1
plot.field.log(gwc_sfc.f09*landfrac.C1, LON.C1, LAT.C1, legend.mar=4, type='def', zlim=c(-2,0.6), col=tim.colors(32)[32:1], Pacific.centric = T)


plot.field.log(gwc_sfc[i.M2,], lon.M2.r, lat.M2, legend.mar=5, type='def', zlim=c(-2.2,0), col=CBR_DRYWET_64[c(3:33,35,37,39,41:62)])
plot.field.log(gwc_sfc[i.M2,], lon.M2.r, lat.M2, legend.mar=4, type='def', zlim=c(-2.2,0), col=CBR_DRYWET_64[c(1,seq(5,64,by=5),64)])
plot.field.log(gwc_sfc[i.M2,], lon.M2.r, lat.M2, legend.mar=5, type='def', zlim=c(-2,-0.2), col=GMT_DRYWET[1:57])



# load CESM lat lon
files.clm = Sys.glob('/Volumes/GoogleDrive/My Drive//CESM_tempfiles/renewmod6_20042005/*.clm*.nc')
nc = nc_open(files.clm[4])#; names(nc$var)
LON.CESM = seq(-180,177.5,by=2.5)
LAT.CESM = ncvar_get(nc, 'lat')
nc_close(nc)

# remember change all NA in the array to 0 before regridding
F_d_LC0.0 = F_d_LC0; F_d_LC0.0[which(is.na(F_d_LC0.0))] = 0
ptm = proc.time()
F_d_LC0.CESM.aw = sp.dissolve(F_d_LC0.0, lon.M2, lat.M2, LON.CESM, LAT.CESM)
proc.time() - ptm
ptm = proc.time()
F_d_LC0.CESM = sp.dissolve.fast(F_d_LC0.0, lon.M2, lat.M2, LON.CESM, LAT.CESM)
proc.time() - ptm

plot.field.log(F_d_LC0.CESM, LON.CESM, LAT.CESM, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
plot.field.log(F_d_LC0.CESM.aw, LON.CESM, LAT.CESM, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

plot.field.log(F_d_LC0.CESM.aw/F_d_LC0.CESM, LON.CESM, LAT.CESM, type='def', zlim=c(-2,2), legend.mar=5, col=TEMP_DIFF_65)


######################

# need to regrid input variables for upscaling
F_d = array(0, dim=dim(C_d))   # using dimension of C_d to define F_d dimension
ind = which(wnd_frc_slt > ustar_it)  # 3-D index for which wind speed > threshold
F_d[ind] = C_tune * C_d[ind] * f_bare[ind] * RHOA[ind] * (wnd_frc_slt[ind]^2-ustar_it[ind]^2)/ustar_it[ind] * (wnd_frc_slt[ind]/ustar_it[ind])^(frag_expt[ind]) * intrmtncy_fct[ind]             # calculate F_d for indices with wind speed > threshold
F_d = sweep(F_d, c(1,2), mss_frc_cly_vld, '*')  # multiply 3-D fields by 2-D fields



##################################################
# plot intermittency 

wt.mean = function(A, B, wt.A, wt.B) {
  if ((dim(A)[1] != dim(B)[1]) | (dim(A)[2] != dim(B)[2])) {stop('dimension of A and B are not equal')}
  sol = array(NaN, dim=dim(A))
  #for (i in 1:(dim(A)[1])) {
  #  for (j in 1:(dim(A)[2])) {
  #    if (is.na(A[i,j])) {sol[i,j] = B[i,j]}
  #    else if (is.na(B[i,j])) {sol[i,j] = A[i,j]}
  #    else {sol[i,j] = (wt.A[i,j]*A[i,j] + wt.B[i,j]*B[i,j])/(wt.A[i,j]+wt.B)[i,j]}
  #{sol[i,j] = (wt.A*A[i,j] + wt.B*B[i,j])/(wt.A+wt.B)}
  #  }
  #}
  
  
  return(sol)
}

num.n0.F_d = function(F_d) {
  num.n0 = array(0, dim=c(dim(F_d)[1], dim(F_d)[2]))
  for (i in 1:(dim(F_d)[1])) {
    num.n0[i,] = rowSums(F_d[i,,] > 0,na.rm=T)
  }
  return(num.n0)
}
num.0.F_d = function(F_d) {
  num.0 = array(0, dim=c(dim(F_d)[1], dim(F_d)[2]))
  for (i in 1:(dim(F_d)[1])) {
    num.0[i,] = rowSums(F_d[i,,] == 0,na.rm=T)
  }
  return(num.0)
}


frag_expt_lim = 3; use.intermittency = T
day=365
# F_d is emission, eta is intermittency, eta_emis is intermittency averaged over time when emission is active
# warning: eta_emis is very big
.env1 = new.env()
#for (d in 1:length(date_vec)) {
#for (d in 23:31) {
for (d in 1:day) {
  print(date_vec[d])
  #filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')
  filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')
  load(filename, envir=.env1)
  ls(.env1)
  if (d==1) {
    F_d = .env1$F_d; eta = .env1$intrmtncy_fct
    F_d.noeta = .env1$F_d.noeta
    eta_emis = F_d / F_d.noeta
    eta_emis[which(F_d==0)]=NaN
    eta_emis = apply(eta_emis, c(1,2), mean, na.rm=T)
    #eta_emis = eta; eta_emis[which(F_d==0)]=NaN; 
    #eta_emis = apply(eta_emis, c(1,2), mean, na.rm=T); 
    #eta_emis[which(is.na(eta_emis))]=0
    #num = num.n0.F_d(F_d)
    
    #eta = apply(.env1$intrmtncy_fct, c(1,2), mean, na.rm=T)
    #plot.field.log(eta, lon.M2, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
    #plot.field.log(eta_emis1, lon.M2, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
    
  }
  else if (d >= 2) {
    F_d1 = .env1$F_d; F_d = F_d + F_d1
    eta = eta + .env1$intrmtncy_fct
    F_d.noeta1 = .env1$F_d.noeta
    F_d.noeta = F_d.noeta + F_d.noeta1
    eta_emis1 = F_d1 / F_d.noeta1
    eta_emis1[which(F_d1==0)]=NaN
    eta_emis1 = apply(eta_emis1, c(1,2), mean, na.rm=T)
    eta_emis = eta_emis + eta_emis1
    #eta_emis1 = .env1$intrmtncy_fct; eta_emis1[which(F_d==0)]=NaN; 
    #eta_emis1 = apply(eta_emis1, c(1,2), mean, na.rm=T); 
    #eta_emis1[which(is.na(eta_emis1))]=0
    #eta_emis = abind(eta_emis, eta_emis1, along=3)
    #eta_emis = eta_emis + eta_emis1
    #eta_emis = wt.mean(eta_emis, eta_emis1, num, num1)
    #num1 = num.n0.F_d(F_d)
    #num = num + num1
  }
}
F_d = apply(F_d, c(1,2), mean, na.rm=T); F_d = F_d/day*365*86400
F_d.noeta = apply(F_d.noeta, c(1,2), mean, na.rm=T); F_d.noeta = F_d.noeta/day*365*86400
eta = apply(eta, c(1,2), mean, na.rm=T); eta = eta/day
eta_emis = eta_emis/day
#eta.F_d = F_d/F_d.noeta
eta_emis = apply(eta_emis, c(1,2), mean, na.rm=T)
#eta_emis = eta_emis/num; eta_emis[which(is.infinite(eta_emis))]=0


plot.field.log(F_d[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR[254:1])
plot.field.log(eta[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-2,0), legend.mar=5, col=WBGYR)
plot.field(eta[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(0,1), col=WBGYR)
plot.field.log(eta_emis[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-2,0), legend.mar=5)
plot.field.log(eta.F_d[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-2,0), legend.mar=5)
plot.field(num[i.M2,], lon.M2.r, lat.M2)


######################################
# current way of plotting intrmittency
day = 365
u = seq(2,12,by=2)
v = 1:6
F_d = eta = array(NaN, dim=c(576, 286, day*length(v))) # warning: 2 large arrays (2.7 GB each)
for (d in 1:day) {
  print(date_vec[d])
  filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')
  
  load(filename, envir=.env1)
  ls(.env1)
  F_d[,,v] = (.env1$F_d)[,,u]
  eta[,,v] = (.env1$intrmtncy_fct)[,,u]
  v = v+6
}



eta_avg = apply(eta, c(1,2), mean, na.rm=T)
plot.field(eta_avg[i.M2,], lon.M2.r, lat.M2, col=WBGYR, type='def', zlim=c(0,1))
plot.field.log(eta_avg[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-2,0), legend.mar=5, col=WBGYR)

eta.dustno0 = eta; eta.dustno0[which(F_d==0)] = NaN
eta.dustno0_avg = apply(eta.dustno0, c(1,2), mean, na.rm=T)
plot.field.log(eta.dustno0_avg[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-2,0), legend.mar=5, col=WBGYR)
plot.field(eta.dustno0_avg[i.M2,], lon.M2.r, lat.M2, col=WBGYR, type='def', zlim=c(0,1))
plot.field.proj(eta.dustno0_avg[i.M2,], lon.M2.r, lat.M2, col=WBGYR)
plot.field.proj(eta.dustno0_avg[i.M2,], lon.M2.r, lat.M2, col=WBGYR, def.zlim=T, zlim=c(0,1))
rm(eta.dustno0) 
rm(eta); rm(F_d)

filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/eta_05x0625_2006_annual_mean.RData', sep='')
save('eta.dustno0_avg','eta_avg',file=filename)

load(filename)
plot.field(eta.dustno0_avg[i.M2,], lon.M2.r, lat.M2, col=WBGYR, type='def', zlim=c(0,1))
plot.field(eta.dustno0_avg[i.M2,], lon.M2.r, lat.M2, col=WBGYR[254:1], type='def', zlim=c(0,1))
A = eta.dustno0_avg; A[which(F_d.L21 <= 1e-4)] = NaN
#plot.field(A[i.M2,], lon.M2.r, lat.M2, col=WBGYR[254:1], type='def', zlim=c(0,1))
plot.field(A[i.M2,], lon.M2.r, lat.M2, col=tim.colors(32)[31:2], type='def', zlim=c(0,0.8))
#plot.field(1-A[i.M2,], lon.M2.r, lat.M2, col=WBGYR, type='def', zlim=c(0,1))
A = eta_avg; A[which(F_d.L21 <= 1e-4)] = NaN
plot.field(A[i.M2,], lon.M2.r, lat.M2, col=tim.colors(32)[31:2], type='def', zlim=c(0,0.8))



##########################################
# plot thresholds from L21
day = 365
#u = seq(6,12,by=6); v = 1:2
.env1 = new.env()
#ustar_ft.wet = ustar_it = array(NaN, dim=c(576, 286, day*length(v))) # warning: 2 large arrays (920 MB each)
for (d in 1:day) {
  print(date_vec[d])
  #filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/thresholds_05x0625_', date_vec[d], '.RData', sep='')
  #filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/thresholds_05x0625_', date_vec[d], '.RData', sep='')
  #save('ustar_ft.wet', 'ustar_it', 'lon.M2', 'lat.M2', file=filename)
  filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006_new/L21_final/thresholds_05x0625_', date_vec[d], '.RData', sep='')
  #save('ustar_it', 'ustar_ft.wet', 'lon.M2', 'lat.M2', file=filename)
  load(filename, envir=.env1)
  ls(.env1)
  #ustar_ft.wet[,,v] = (.env1$ustar_ft.wet)[,,u]
  #ustar_it[,,v] = (.env1$ustar_it)[,,u]
  #v = v+6
  if (d==1) {
    ustar_ft.wet = .env1$ustar_ft.wet
    ustar_it = .env1$ustar_it
  }
  else if (d >= 2) {
    ustar_ft.wet = ustar_ft.wet + .env1$ustar_ft.wet
    ustar_it = ustar_it + .env1$ustar_it
  }
}
ustar_it = apply(ustar_it, c(1,2), mean, na.rm=T) / day
ustar_ft.wet = apply(ustar_ft.wet, c(1,2), mean, na.rm=T) / day


plot.field(ustar_it[i.M2,], lon.M2.r, lat.M2)
plot.field(ustar_ft.wet[i.M2,], lon.M2.r, lat.M2)
plot.field((ustar_ft.wet/ustar_it)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(1,4.5))
plot.field((ustar_ft.wet/ustar_it)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(1,4.5), col=CBR_DRYWET_65)
plot.field((ustar_ft.wet/ustar_it)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(1,4.5), col=CBR_DRYWET_65[c(4:61)])


plot.field.proj(ustar_it[i.M2,], lon.M2.r, lat.M2)

directory = 'L21_final'

filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006_new/', directory, '/thresholds_05x0625_2006_annual_mean.RData', sep='')
save('ustar_it', 'ustar_ft.wet', 'lon.M2', 'lon.M2.r', 'lat.M2', file=filename) 
#filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/thresholds_05x0625_2006_annual_mean.RData', sep='')
#save('ustar_it','ustar_ft.wet',file=filename)


######################################################################
# a quick comparison between K21 DUSTCOMM emissions and our emissions

library(ggplot2)

nc = nc_open('DustCOMM_dust_emission_flux_annual_PM20.nc')
names(nc$var)
F_d.K21 = ncvar_get(nc, 'Mean')
lon.K21 = ncvar_get(nc, 'lon')
lat.K21 = ncvar_get(nc, 'lat')
area.K21 = get.area(lon.K21, lat.K21, 'm2')
sum(F_d.K21*area.K21, na.rm=TRUE)/1e9   # Tg/yr-1
#plot.field.log(F_d.K21/5693*5000, lon.K21, lat.K21, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
indj.K = 18:93; indi.K = c(5:144,1:4)
lon.K21.r = lon.K21[indi.K]; lon.K21.r[141:144]=lon.K21.r[141:144]+360
lat.K21.r = lat.K21[indj.K]
#plot.field.log(F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
plot.field.log(F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR)


# Jasper says remove the Northern US emissions (7 July 2022)
F_d.K21[which(lon.K21<(-30)), which(lat.K21>45)]=0
range(F_d.K21[which(lon.K21<(-30)), which(lat.K21>50)])
plot.field.log(F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR); frame.region(frame='source',lwd=1.5)


# extract K21a region data
K21.region = extract.value.region(F_d.K21[indi.K,indj.K]/5679*5000, lon.K21.r, lat.K21.r, output.matrix.3d=TRUE, frame='source')
K21.region.grid = K21.region$matrix.3d  # gridded emissions separated by regions
plot.field.log(K21.region.grid[,,1], lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR)

area.K21 = get.area(lon.K21.r, lat.K21.r, 'm2')
F_d.K21.reg = apply(sweep(K21.region.grid, c(1,2), area.K21/1e9, '*'), c(3), sum, na.rm=T)

# change to Bullard et al. (2016) for high-latitude dust emission
F_d.K21.reg[10] = 90*5/2


# also K21b Table 1 regional dust
print(name.region.Ridley.source.K21)
#F_d.K21b.reg = c(0.88,0.72,0.56,1.38,0.58,0.13,0.16,0.19,0.10)*1e3*(47/50)
#F_d.K21b.reg.negsig = c(0.64,0.47,0.15,0.97,0.42,0.05,0.09,0.10,0.06)*1e3*(47/50)
#F_d.K21b.reg.possig = c(1.44,1.11,1.77,2.59,1.12,0.24,0.29,0.35,0.19)*1e3*(47/50)
F_d.K21b.reg = c(18,16,13,29,13,3,3,4,2)*(5000/100)
F_d.K21b.reg.negsig = c(14,7,4,27,10,1,2,2,1)*(5000/100)
F_d.K21b.reg.possig = c(22,21,20,32,15,4,5,6,3)*(5000/100)

# add Bullard et al. (2016) high latitude dust emissions
F_d.K21b.reg = c(F_d.K21b.reg, 90*5/2)
F_d.K21b.reg.negsig = c(F_d.K21b.reg.negsig,0)
F_d.K21b.reg.possig = c(F_d.K21b.reg.possig,0)

F_d.K21b.reg.negsd = F_d.K21b.reg.negsig - F_d.K21b.reg
F_d.K21b.reg.possd = F_d.K21b.reg.possig - F_d.K21b.reg

#F_d.K21.05 = sp.smooth(F_d.K21[indi.K,indj.K], lon.K21.r, lat.K21.r, lon.M2.r, lat.M2)
#sum(F_d.K21.05*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
#plot.field.log(F_d.K21.05/5679*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
#plot.field.log(F_d.L21[i.M2,]/8778*5000 / F_d.K21.05, lon.M2.r, lat.M2, type='def', zlim=c(-1,1), legend.mar=5, col=TEMP_DIFF_65)

# par(mai=c(0.6,0.6,0.1,0.1))
# plot(log10(F_d.L21[i.M2,]/8778*5000), log10(F_d.K21.05/5679*5000), 'p', pch=17, col='blue', xlim=c(-7,1), ylim=c(-7,1))
# A = mat_2dto1d(F_d.L21[i.M2,]/8778*5000)
# B = mat_2dto1d(F_d.K21.05/5679*5000)
# cor(A, B, 'complete.obs')
# 
# plot(log10(F_d.K14_def[i.M2,]/25187*5000), log10(F_d.K21.05/5679*5000), 'p', pch=17, col='blue', xlim=c(-7,1), ylim=c(-7,1))
# C = mat_2dto1d(F_d.K14_def[i.M2,]/25187*5000)
# D = mat_2dto1d(F_d.K21.05/5679*5000)
# cor(C, D, 'complete.obs')

# regrid and extract region
F_d.L21.re = sp.dissolve(F_d.L21[i.M2,], lon.M2.r, lat.M2, lon.K21.r, lat.K21.r)
plot.field.log(F_d.L21.re/11748*5000, lon.K21.r, lat.K21.r, legend.mar=5)
# Patagonia emissions
F_d.L21.Pat = extract.value.selfdef.region(F_d.L21.re/11748*5000, lon.K21.r, lat.K21.r, low.lat=(-58), high.lat=(-38), left.lon=(-80), right.lon=(-65), output.matrix=TRUE)
F_d.L21.Pat = sum(F_d.L21.Pat$matrix.2d * area.K21/1e9, na.rm=T)  # total Patagonia emissions in Tg / yr
#plot.field.log(F_d.L21.Pat$matrix.2d, lon.K21.r, lat.K21.r, legend.mar=5)
L21.region = extract.value.region(F_d.L21.re/11748*5000, lon.K21.r, lat.K21.r, output.matrix.3d=TRUE, frame='source')
names(L21.region)
L21.region.grid = L21.region$matrix.3d  # gridded emissions separated by regions
F_d.L21.reg = apply(sweep(L21.region.grid, c(1,2), area.K21/1e9, '*'), c(3), sum, na.rm=T)
F_d.L21.reg[10] = F_d.L21.reg[10] + F_d.L21.Pat


F_d.K14_def.re = sp.dissolve(F_d.K14_def[i.M2,], lon.M2.r, lat.M2, lon.K21.r, lat.K21.r)
plot.field.log(F_d.K14_def.re/29254*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR)
K14.region = extract.value.region(F_d.K14_def.re/29254*5000, lon.K21.r, lat.K21.r, output.matrix.3d=TRUE, frame='source')
K14.region.grid = K14.region$matrix.3d  # gridded emissions separated by regions
plot.field.log(K14.region.grid[,,1], lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR)
F_d.K14.reg = apply(sweep(K14.region.grid, c(1,2), area.K21/1e9, '*'), c(3), sum, na.rm=T)
# Patagonia emissions
F_d.K14.Pat = extract.value.selfdef.region(F_d.K14_def.re/29254*5000, lon.K21.r, lat.K21.r, low.lat=(-58), high.lat=(-39), left.lon=(-80), right.lon=(-65), output.matrix=TRUE)
F_d.K14.Pat = sum(F_d.K14.Pat$matrix.2d * area.K21/1e9, na.rm=T)  # total Patagonia emissions in Tg / yr
F_d.K14.reg[10] = F_d.K14.reg[10] + F_d.K14.Pat



# comparing L21 against K21
plot.field.log(F_d.L21[i.M2,]/11748*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

plot.field.log(F_d.L21.re/11748*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR)

par(mai=c(0.6,0.6,0.1,0.1))
plot(log10(F_d.L21.re/11748*5000), log10(F_d.K21[indi.K,indj.K]/5679*5000), 'p', pch=17, col='blue', xlim=c(-5,0.5), ylim=c(-5,0.5))
A = mat_2dto1d(log10(F_d.L21.re/11748*5000)); A[which(is.infinite(A))] = NaN; A[which(A< (-5))] = NaN
B = mat_2dto1d(log10(F_d.K21[indi.K,indj.K]/5679*5000)); B[which(is.infinite(B))] = NaN; B[which(B< (-5))] = NaN
cor(A, B, 'complete.obs')


plot.field.log(F_d.K21[indi.K,indj.K]/5679*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR); frame.region(frame='source',lwd=1.5)
CC = F_d.L21.re/11748*5000 / (F_d.K21[indi.K,indj.K]/5679*5000); CC[which(is.infinite(CC))] = NaN; CC[which(F_d.K21[indi.K,indj.K]/5679*5000 < 1e-4)] = NaN
plot.field.log(CC, lon.K21.r, lat.K21.r, type='def', zlim=c(-1,1), legend.mar=3, col=TEMP_DIFF_65)



par(mai=c(0.6,0.6,0.1,0.1))
plot(log10(F_d.L21.re/11748*5000), log10(F_d.K21[indi.K,indj.K]/5679*5000), 'p', pch=17, col='blue', xlim=c(-4,0.3), ylim=c(-4,0.3))
A = mat_2dto1d(log10(F_d.L21.re/11748*5000)); A[which(is.infinite(A))] = NaN; A[which(A< (-4))] = NaN
B = mat_2dto1d(log10(F_d.K21[indi.K,indj.K]/5679*5000)); B[which(is.infinite(B))] = NaN; B[which(B< (-4))] = NaN
cor(A, B, 'complete.obs'); sqrt(mean((A - B)^2, na.rm=T))

plot.field.log(L21.region.grid[,,1], lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR)


cor(mat_2dto1d(F_d.L21.re/11748*5000), mat_2dto1d(F_d.K21[indi.K,indj.K]/5679*5000), 'complete.obs')
cor(mat_2dto1d(F_d.K14_def.re/29254*5000), mat_2dto1d(F_d.K21[indi.K,indj.K]/5679*5000), 'complete.obs')
cor(mat_2dto1d(F_d.Z03.Z.01.re/522), mat_2dto1d(F_d.K21[indi.K,indj.K]/5679*5000), 'complete.obs')
cor(mat_2dto1d(F_d.Z03.G.re/442), mat_2dto1d(F_d.K21[indi.K,indj.K]/5679*5000), 'complete.obs')

# plot gridded data by region
par(mai=c(0.6,0.6,0.1,0.1), ps=16, mgp=c(1.8,0.5,0))
plot(NA,NA, xlim=c(-2,0), ylim=c(-2,0), xaxt='n', yaxt='n', xlab="this study's emissions", ylab='K21a emissions')
axis(side = 1, at=((-4):0), labels=expression('10'^-4,'10'^-3,'10'^-2,'10'^-1,'0')); axis(side = 2, at=((-4):0), labels=expression('10'^-4,'10'^-3,'10'^-2,'10'^-1,'0'))
abline(a=0,b=1, lwd=2); abline(a=-1,b=1, lwd=2, lty=2); abline(a=-2,b=1, lwd=2, lty=3); abline(a=1,b=1, lwd=2, lty=2); abline(a=2, b=1, lwd=2, lty=3)
abline(v=-4, col='grey', lwd=2); abline(h=-4, col='grey', lwd=2)
.col = c('#8AD2D8','#C6A68E','#558AA4','#F15E3D','#56704B','#CC3B7C','#005594','#89B65A','#EBC98E','#EFCB5E')
for (n in 2) {
  #points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=n+1, pch=19); points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=1)
  #points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=n+1, pch=8)
  #points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=.col, pch=19); points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=1, lwd=1.5)
  #points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=.col, pch=8)
  points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=.col[6], pch=19, cex=1.3); points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=1, cex=1.3, lwd=1.5)
  
}

legend('center', col=.col, legend=c(name.region.Ridley.source.K21, 'RestofWorld'), pch=19)

n=7
cor(mat_2dto1d(L21.region.grid[,,n]) , mat_2dto1d(K21.region.grid[,,n]) )

cor(mat_2dto1d(K14.region.grid[,,n]) , mat_2dto1d(K21.region.grid[,,n]) )

cor(mat_2dto1d(Z03.G.region.grid[,,n]) , mat_2dto1d(K21.region.grid[,,n]) )
cor(mat_2dto1d(Z03.Z.region.grid[,,n]) , mat_2dto1d(K21.region.grid[,,n]) )


# regional scatterplot

par(mai=c(0.6,0.6,0.1,0.1), mgp=c(1.6,0.4,0), ps=13)
plot(NA,NA, xlim=c(40,1580), ylim=c(40,1580), xlab="K14 (red) and our study's (blue) emissions (Tg/yr)", ylab='K21 emissions (Tg/yr)')
abline(a=0,b=1,lwd=2); abline(v=400,lwd=1,col='grey'); abline(v=800,lwd=1,col='grey'); abline(v=1200,lwd=1,col='grey'); abline(h=400,lwd=1,col='grey'); abline(h=800,lwd=1,col='grey'); abline(h=1200,lwd=1,col='grey')

points(F_d.L21.reg, F_d.K21.reg, pch=c(1:10), 'p', col='blue', cex=2.5, lwd=2)
points(F_d.K14.reg, F_d.K21.reg, pch=c(1:10), 'p', col='red', cex=2.5, lwd=2)
arrows(F_d.L21.reg, F_d.K21.reg+F_d.K21b.reg.negsd, F_d.L21.reg, F_d.K21.reg+F_d.K21b.reg.possd, length=0.1, angle=90, code=3, col='blue', lwd=2)
arrows(F_d.K14.reg, F_d.K21.reg+F_d.K21b.reg.negsd, F_d.K14.reg, F_d.K21.reg+F_d.K21b.reg.possd, length=0.1, angle=90, code=3, col='red', lwd=2)

cor(F_d.L21.reg, F_d.K21.reg)^2; rmse(F_d.L21.reg, F_d.K21.reg)
cor(F_d.K14.reg, F_d.K21.reg)^2; rmse(F_d.K14.reg, F_d.K21.reg)

points(F_d.L21.reg, F_d.K21b.reg, pch=c(1:10), 'p', col='blue', cex=2.5, lwd=2)
points(F_d.K14.reg, F_d.K21b.reg, pch=c(1:10), 'p', col='red', cex=2.5, lwd=2)
arrows(F_d.L21.reg, F_d.K21b.reg+F_d.K21b.reg.negsd, F_d.L21.reg, F_d.K21b.reg+F_d.K21b.reg.possd, length=0.1, angle=90, code=3, col='blue', lwd=2)
arrows(F_d.K14.reg, F_d.K21b.reg+F_d.K21b.reg.negsd, F_d.K14.reg, F_d.K21b.reg+F_d.K21b.reg.possd, length=0.1, angle=90, code=3, col='red', lwd=2)

cor(F_d.L21.reg, F_d.K21b.reg)^2; rmse(F_d.L21.reg, F_d.K21b.reg)
cor(F_d.K14.reg, F_d.K21b.reg)^2; rmse(F_d.K14.reg, F_d.K21b.reg)


#legend('center',legend=c(name.region.Ridley.source.K21,'RestofWorld'), pch=1:10)
legend('center',legend=c(name.region.Ridley.source.K21,'high-latitude'), pch=1:10)

##########
# regional scatterplot in log scale

par(mai=c(0.6,0.6,0.1,0.1), mgp=c(1.6,0.4,0), ps=13)
plot(NA,NA, xlim=c(20,3000), ylim=c(20,3000), xlab="K14 (red) and our study's (blue) emissions (Tg/yr)", ylab='K21 emissions (Tg/yr)', log='xy')
abline(a=0,b=1,lwd=2); abline(v=10,lwd=1,col='grey'); abline(v=100,lwd=1,col='grey'); abline(v=1000,lwd=1,col='grey'); abline(h=10,lwd=1,col='grey'); abline(h=100,lwd=1,col='grey'); abline(h=1000,lwd=1,col='grey')

points(F_d.L21.reg, F_d.K21b.reg, pch=c(1:10), 'p', col='blue', cex=2.5, lwd=2)
points(F_d.K14.reg, F_d.K21b.reg, pch=c(1:10), 'p', col='red', cex=2.5, lwd=2)
arrows(F_d.L21.reg, F_d.K21b.reg+F_d.K21b.reg.negsd, F_d.L21.reg, F_d.K21b.reg+F_d.K21b.reg.possd, length=0.1, angle=90, code=3, col='blue', lwd=2) # bug found 12 Jan 2023: it should be F_d.K21b.reg+F_d.K21b.reg.possd instead of F_d.K21.reg+F_d.K21b.reg.possd
arrows(F_d.K14.reg, F_d.K21b.reg+F_d.K21b.reg.negsd, F_d.K14.reg, F_d.K21b.reg+F_d.K21b.reg.possd, length=0.1, angle=90, code=3, col='red', lwd=2) # bug found 12 Jan 2023: it should be F_d.K21b.reg+F_d.K21b.reg.possd instead of F_d.K21.reg+F_d.K21b.reg.possd

cor(F_d.L21.reg, F_d.K21b.reg)^2; rmse(F_d.L21.reg, F_d.K21b.reg)
cor(F_d.K14.reg, F_d.K21b.reg)^2; rmse(F_d.K14.reg, F_d.K21b.reg)


###
# comparing K14 against K21
plot.field.log(F_d.K14_def[i.M2,]/29254*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

plot.field.log(F_d.K14_def.re/29254*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR)
plot.field.log(F_d.K14_def.re/29254*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR, horizontal = T)
plot.field.log(F_d.K21[indi.K,indj.K]/5679*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5)
CC = F_d.K14_def.re/29254*5000 / (F_d.K21[indi.K,indj.K]/29254*5000); CC[which(is.infinite(CC))] = NaN; CC[which(F_d.K21[indi.K,indj.K]/29254*5000 < 1e-4)] = NaN
plot.field.log(CC, lon.K21.r, lat.K21.r, type='def', zlim=c(-1,1), legend.mar=3, col=TEMP_DIFF_65)


par(mai=c(0.6,0.6,0.1,0.1))
plot(log10(F_d.K14_def.re/29254*5000), log10(F_d.K21[indi.K,indj.K]/5679*5000), 'p', pch=17, col='blue', xlim=c(-4,0.3), ylim=c(-4,0.3))
A = mat_2dto1d(log10(F_d.K14_def.re/29254*5000)); A[which(is.infinite(A))] = NaN; A[which(A< (-4))] = NaN
B = mat_2dto1d(log10(F_d.K21[indi.K,indj.K]/5679*5000)); B[which(is.infinite(B))] = NaN; B[which(B< (-4))] = NaN
cor(A, B, 'complete.obs'); sqrt(mean((A - B)^2, na.rm=T))



# plot by region
par(mai=c(0.6,0.6,0.1,0.1), ps=16, mgp=c(1.8,0.5,0))
plot(NA,NA, xlim=c(-2,0), ylim=c(-2,0), xaxt='n', yaxt='n', xlab="K14 emissions", ylab='K21 emissions')
axis(side = 1, at=((-4):0), labels=expression('10'^-4,'10'^-3,'10'^-2,'10'^-1,'0')); axis(side = 2, at=((-4):0), labels=expression('10'^-4,'10'^-3,'10'^-2,'10'^-1,'0'))
abline(a=0,b=1, lwd=2); abline(a=-1,b=1, lwd=2, lty=2); abline(a=-2,b=1, lwd=2, lty=3); abline(a=1,b=1, lwd=2, lty=2); abline(a=2, b=1, lwd=2, lty=3)
abline(v=-4, col='grey', lwd=2); abline(h=-4, col='grey', lwd=2)
.col = c('#8AD2D8','#C6A68E','#558AA4','#F15E3D','#56704B','#CC3B7C','#005594','#89B65A','#EBC98E','#EFCB5E')
for (n in 1:10) {
  #points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=n+1, pch=19); points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=1)
  #points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=n+1, pch=8)
  #points(log10(K14.region.grid[,,n]), log10(K21.region.grid[,,n]), col=.col, pch=19); points(log10(K14.region.grid[,,n]), log10(K21.region.grid[,,n]), col=1, lwd=1.5)
  points(log10(K14.region.grid[,,n]), log10(K21.region.grid[,,n]), col=.col[6], pch=19, cex=1.3); points(log10(K14.region.grid[,,n]), log10(K21.region.grid[,,n]), col=1, cex=1.3, lwd=1.5)
  #points(log10(L21.region.grid[,,n]), log10(K21.region.grid[,,n]), col=.col, pch=8)
}



#######################################################
# get MERRA2 time series for a grid
# function for getting MERRA2 time series
get.merra2.ts = function(files, date_vec, dt, latlonbound=rbind(c(16,18),c(17,19)), var=NULL) {
  # get annual dust emission total per grid (kg m-2 yr-1) 
  strings = strsplit(files, '[.]')
  date_files = as.double(lapply(strings, `[[`, 3)) # date
  
  # get indices within the range of lat long
  indx = which(lon.M2 >= latlonbound[1,1] & lon.M2 <= latlonbound[1,2])
  indy = which(lat.M2 >= latlonbound[2,1] & lat.M2 <= latlonbound[2,2])
  
  ind = match(date_vec, date_files)
  print(ind)
  .env1 = new.env()
  VAR = NULL
  for (d in 1:length(date_vec)) {
    print(date_vec[d])
    nc = nc_open(files[ind[d]])#; names(nc$var)#, start=c(11,26), count=c(5,5)
    
    #if (d==1) {VAR = ncvar_get(nc, var, start=c(i1,j1,t1), count=c(il,jl,tl))} # m/s
    #else if (d >= 2) {VAR = VAR + ncvar_get(nc, var, start=c(i1,j1,t1), count=c(il,jl,tl))}
    VAR = c(VAR, ncvar_get(nc, var, start=c(indx[1],indy[1],1), count=c((indx[length(indx)]),(indy[length(indy)]),24))[,,indt]) # m/s
    
  }
  #dim(VAR)
  
  return(list(VAR=VAR))
}

date_vec = make.date.vec(20060101,20060331); directory = 'L21_final'
ustar_ft.wet = ustar_it = PBLH = Obu_L = wnd_frc_slt = intrmtncy_fct = NULL
.env1 = new.env()
for (d in 1:length(date_vec)) {
  print(date_vec[d])
  #filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006_new/', directory, '/intermittency_05x0625_', date_vec[d], '.RData', sep='')
  filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006_new/', directory, '/intermittency_05x0625_', date_vec[d], '.RData', sep='')
  
  load(filename, envir=.env1)
  Obu_L = c(Obu_L, apply((.env1$Obu_L)[indx,indy,], 3, mean, na.rm=T))
  PBLH = c(PBLH, apply((.env1$PBLH)[indx,indy,], 3, mean, na.rm=T))
  intrmtncy_fct = c(intrmtncy_fct, apply((.env1$intrmtncy_fct)[indx,indy,], 3, mean, na.rm=T))
  wnd_frc_slt = c(wnd_frc_slt, apply((.env1$wnd_frc_slt)[indx,indy,], 3, mean, na.rm=T))
  ustar_it = c(ustar_it, apply((.env1$ustar_it)[indx,indy,], 3, mean, na.rm=T))
  ustar_ft.wet = c(ustar_ft.wet, apply((.env1$ustar_ft.wet)[indx,indy,], 3, mean, na.rm=T))
  #c('Obu_L', 'PBLH', 'intrmtncy_fct', 'ustar_ft.wet', 'ustar_it', 'wnd_frc_slt', 'lon.M2', 'lat.M2')
}

par(mai=c(0.6,0.6,0.1,0.1))
plot(wnd_frc_slt, intrmtncy_fct, 'p', col='blue', pch=19, xlim=c(0.016,0.52))
abline(v=mean(ustar_it), lwd=2, lty=2)
abline(v=mean(ustar_ft.wet), lwd=2, lty=2)


plot(PBLH, intrmtncy_fct, 'p', col='blue', pch=19)
plot(PBLH, wnd_frc_slt, 'p', col='blue', pch=19)
library(ggplot2)


limm = 100; k = 0.387
zeta = (-PBLH/Obu_L); zeta[which(zeta>limm)] = limm; zeta[which(zeta<(-limm))] = -limm
sigma = wnd_frc_slt * (12 - 0.5 * PBLH/Obu_L)^0.333; sigma[which(is.na(sigma))] = 0
u01 = wnd_frc_slt / k * log(0.1/1e-4)
uit01 = ustar_it / k * log(0.1/1e-4)
uft01 = ustar_ft.wet / k * log(0.1/1e-4)
intrmtncy_fct1 = intrmtncy_fct; intrmtncy_fct1[which(intrmtncy_fct>0.7 & wnd_frc_slt<0.2)] = NaN
dataframe = data.frame(PBLH=PBLH, wnd_frc_slt=wnd_frc_slt, intrmtncy_fct=intrmtncy_fct1, zeta=zeta, sigma=sigma, u01=u01, uft01=uft01, uit01=uit01)
ggplot(dataframe, aes(x = u01, y = intrmtncy_fct, color = zeta)) + geom_point(size = 2.5) + scale_colour_gradientn(colours = TEMP_DIFF_19[c(3:8,12:17)], limits=c(-limm, limm)) + xlab("hourly mean wind speed (u, m/s) at 0.1 m") + ylab("hourly mean intermittency") + geom_vline(xintercept=mean(uft01), lty=2, lwd=1.1) + geom_vline(xintercept=mean(uit01), lty=2, lwd=1.1) + xlim(0,10)#+ theme(panel.background = element_rect(fill = "lightgrey",  colour = "lightgrey", size = 0.5, linetype = "solid"))

library(RColorBrewer)
Spe = brewer.pal(n = 11, name = "Spectral")
#Spe = Spe[length(RdBu):1]
getPalette = colorRampPalette(Spe)
getPalette(32)

ggplot(dataframe, aes(x = u01, y = intrmtncy_fct, color = zeta)) + geom_point(size = 2.5) + scale_colour_gradientn(colours = getPalette(16), limits=c(-limm, limm)) + xlab("hourly mean wind speed (u, m/s) at 0.1 m") + ylab("hourly mean intermittency") + geom_vline(xintercept=mean(uft01), lty=2, lwd=1.1) + geom_vline(xintercept=mean(uit01), lty=2, lwd=1.1) + xlim(0,10)#+ theme(panel.background = element_rect(fill = "white",  colour = "lightblue", size = 0.5, linetype = "solid"))


ggplot(dataframe, aes(x = u01, y = intrmtncy_fct, color = zeta)) + geom_point(size = 2.5) + scale_colour_gradientn(fill = TEMP_DIFF_19[c(3:9,11:17)], colors="black", pch=21, limits=c(-limm, limm)) + xlab("hourly mean wind speed (u, m/s) at 0.1 m") + ylab("hourly mean intermittency") + geom_vline(xintercept=mean(uft01), lty=2, lwd=1.1) + geom_vline(xintercept=mean(uit01), lty=2, lwd=1.1) + xlim(0,10) 

ggplot(dataframe, aes(x = u01, y = intrmtncy_fct, color = sigma)) + geom_point(size = 2.5) + scale_colour_gradientn(colours = tim.colors(32), limits=c(0,1.5)) + xlab("hourly mean wind speed (u, m/s) at 0.1 m") + ylab("hourly mean intermittency") + geom_vline(xintercept=mean(uft01), lty=2, lwd=1.1) + geom_vline(xintercept=mean(uit01), lty=2, lwd=1.1) + xlim(0,10) 



plot(1:200, PBLH[1:200], 'l')
lines(1:200, 500*intrmtncy_fct[1:200], 'l', col='blue')
lines(1:200, 500*intrmtncy_fct[1:200], 'l', col='blue')

zeta=(-PBLH/Obu_L)


.col=tim.colors(32)

range(PBLH)
seq(round(63.08811),round(3193.62183),length=length(.col))

color.var <- cut(PBLH, breaks=seq(round(63.08811),round(3193.62183),length=(length(.col)+1)), labels=.col)


breakk=seq(round(63.08811),round(3193.62183),length=(length(.col)+1))
for (i in 1:length(.col)) {
  data$Colour[data$col_name2>=3]="red"
}

par(mai=c(0.6,0.6,0.1,0.1))
plot(wnd_frc_slt, intrmtncy_fct, 'p', col=color.var, pch=19, xlim=c(0.016,0.52))
abline(v=mean(ustar_it), lwd=2, lty=2)
abline(v=mean(ustar_ft.wet), lwd=2, lty=2)

