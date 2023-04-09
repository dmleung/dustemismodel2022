########################################################
########################################################
# 12 Jan 2023
# script for response to reviewers' comments
########################################################
########################################################

setwd('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR')
# load libraries 
library(ncdf4); library(fields); library(maps); library(Metrics); library(pracma); library(abind)
# load functions
source("get_geo.R"); source("get_met.R"); source("sptial_plot_fns.R")
source("dust_research_fns.R")
# load color scheme
load("WBGYR_scheme.RData"); load('TEMP_DIFF.RData'); load('LEAF_scheme.RData'); load('CBR_DRYWET.RData'); load('BAGYOR_scheme.RData')

# indices
indi.M2 = 1:576; indj.M2 = 65:350 
indi.M05 = 1:540; indj.M05 = 65:350 #not sure why MERRA2 grid is not available. 12 Jan 2023

########################################################
# MERRA-2 dust emissions

#files = Sys.glob('/Users/dannymleung/Google Drive/My Drive/MERRA2/MERRA-2_dust_var/*')
files = Sys.glob('/Volumes/GoogleDrive/My Drive/MERRA2/MERRA-2_dust_var/*')

# get long lat
nc = nc_open(files[1])
lon.M05 = ncvar_get(nc, 'lon'); lat.M05 = ncvar_get(nc, 'lat')
length(lon.M05); length(lat.M05)
lon.M05 = ncvar_get(nc, 'lon')[indi.M05]; lat.M05 = ncvar_get(nc, 'lat')[indj.M05]
nc_close(nc)
i.M05 = c(17:540,1:16)
lon.M05.r = lon.M05[i.M05]; lon.M05.r[525:540]=lon.M05.r[525:540]+360

# get area 2-D map
area.M05 = get.area(lon.M05.r, lat.M05, 'm2')

# get dust emissions
DUEM = array(NaN, dim=c(length(lon.M05), length(lat.M05), length(files)))
for (m in 1:length(files)) {
  print(m)
  nc = nc_open(files[m])
  names(nc$var)
  DUEM[,,m] = (ncvar_get(nc, 'DUEM001') + ncvar_get(nc, 'DUEM002') + ncvar_get(nc, 'DUEM003') + ncvar_get(nc, 'DUEM004') + ncvar_get(nc, 'DUEM005'))[indi.M05,indj.M05]  # kg m-2 s-1
}

# plot emissions
DUEM_annavg = apply(DUEM, c(1,2), mean, na.rm=T)*86400*365
sum(DUEM_annavg*area.M05, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(DUEM_annavg[i.M05,]*5000/1561, lon.M05.r, lat.M05, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)#; frame.region(frame='source', lwd=1.5)
plot.field.log(DUEM_annavg[i.M05,]*5000/1561, lon.M05.r, lat.M05, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR, horizontal=T, legend.only=T)

# annual mean LAI (use get.merra2.field from MERRA2_dustmod_05x06_analysis_new.R)
FF = get.merra2.field(files_lnd, date_vec, dt, var='LAI')
LAI.ann = FF$VAR; dim(LAI.ann)
LAI.filt.arid = LAI.ann; LAI.filt.arid[which(LAI.filt.arid>1)] = NaN; LAI.filt.arid[which(LAI.filt.arid <= 1)] = 1
LAI.filt.nonarid = LAI.ann; LAI.filt.nonarid[which(LAI.filt.nonarid <= 1)] = NaN; LAI.filt.nonarid[which(LAI.filt.nonarid > 1)] = 1

DUEM.M2 = sp.dissolve(DUEM_annavg, lon.M05.r, lat.M05, lon.M2.r, lat.M2)
sum(DUEM.M2*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
sum(DUEM.M2*LAI.filt.nonarid*area.M2, na.rm=TRUE)/1e9   # nonarid Tg/yr-1
sum(DUEM.M2*LAI.filt.arid*area.M2, na.rm=TRUE)/1e9   # arid Tg/yr-1


# compare with DustCOMM (K21) emissions
nc = nc_open('/Volumes/GoogleDrive/My Drive/CESM/CESMDUSTinR/DustCOMM_dust_emission_flux_annual_PM20.nc'); names(nc$var)
F_d.K21 = ncvar_get(nc, 'Mean')
lon.K21 = ncvar_get(nc, 'lon'); lat.K21 = ncvar_get(nc, 'lat')
area.K21 = get.area(lon.K21, lat.K21, 'm2')
sum(F_d.K21*area.K21, na.rm=TRUE)/1e9   # Tg/yr-1
indj.K = 18:93; indi.K = c(5:144,1:4)
lon.K21.r = lon.K21[indi.K]; lon.K21.r[141:144]=lon.K21.r[141:144]+360; lat.K21.r = lat.K21[indj.K]
plot.field.log(F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
# Jasper says remove the Northern US emissions (7 July 2022)
F_d.K21[which(lon.K21<(-30)), which(lat.K21>45)]=0
range(F_d.K21[which(lon.K21<(-30)), which(lat.K21>50)])
plot.field.log(F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR)

# use K21b Table for regional emissions
name.region.Ridley.source
F_d.K21b.reg = c(18,16,13,29,13,3,3,4,2)*(5000/100)
F_d.K21b.reg.negsig = c(14,7,4,27,10,1,2,2,1)*(5000/100)
F_d.K21b.reg.possig = c(22,21,20,32,15,4,5,6,3)*(5000/100)
# add Bullard et al. (2016) high latitude dust emissions
F_d.K21b.reg = c(F_d.K21b.reg, 90*5/2)
F_d.K21b.reg.negsig = c(F_d.K21b.reg.negsig,0)
F_d.K21b.reg.possig = c(F_d.K21b.reg.possig,0)
F_d.K21b.reg.negsd = F_d.K21b.reg.negsig - F_d.K21b.reg
F_d.K21b.reg.possd = F_d.K21b.reg.possig - F_d.K21b.reg



# regrid MERRA-2 emissions to DustCOMM (K21) and extract region 
# code got from MERRA2_dustmod_05_06_analysis_new.R. 12 Jan 2023
DUEM.re = sp.dissolve(DUEM_annavg[i.M05,], lon.M05.r, lat.M05, lon.K21.r, lat.K21.r)
area.K21 = get.area(lon.K21.r, lat.K21.r, 'm2')
sum(DUEM.re*area.K21, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(DUEM.re/1567*5000, lon.K21.r,lat.K21.r, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
# Patagonia emissions
DUEM.Pat = extract.value.selfdef.region(DUEM.re/1567*5000, lon.K21.r, lat.K21.r, low.lat=(-58), high.lat=(-38), left.lon=(-80), right.lon=(-65), output.matrix=TRUE)
DUEM.Pat = sum(DUEM.re.Pat$matrix.2d * area.K21/1e9, na.rm=T)  # total Patagonia emissions in Tg / yr
#plot.field.log(DUEM.Pat$matrix.2d, lon.K21.r, lat.K21.r, legend.mar=5)
DUEM.region = extract.value.region(DUEM.re/1567*5000, lon.K21.r, lat.K21.r, output.matrix.3d=TRUE, frame='source')
names(DUEM.region)
DUEM.region.grid = DUEM.region$matrix.3d  # gridded emissions separated by regions
DUEM.reg = apply(sweep(DUEM.region.grid, c(1,2), area.K21/1e9, '*'), c(3), sum, na.rm=T)
DUEM.reg[10] = DUEM.reg[10] + DUEM.Pat

# regional scatterplot in log scale

par(mai=c(0.6,0.6,0.1,0.1), mgp=c(1.6,0.4,0), ps=13)
#plot(NA,NA, xlim=c(20,3000), ylim=c(20,3000), xlab="K14 (red) and our study's (blue) emissions (Tg/yr)", ylab='K21 emissions (Tg/yr)', log='xy')
plot(NA,NA, xlim=c(20,3000), ylim=c(20,3000), xlab="MERRA-2 emissions (Tg/yr)", ylab='DustCOMM emissions (Tg/yr)', log='xy')
abline(a=0,b=1,lwd=2); abline(v=10,lwd=1,col='grey'); abline(v=100,lwd=1,col='grey'); abline(v=1000,lwd=1,col='grey'); abline(h=10,lwd=1,col='grey'); abline(h=100,lwd=1,col='grey'); abline(h=1000,lwd=1,col='grey')

points(DUEM.reg, F_d.K21b.reg, pch=c(1:10), 'p', col='darkorange4', cex=2.5, lwd=2)
#points(F_d.K14.reg, F_d.K21b.reg, pch=c(1:10), 'p', col='red', cex=2.5, lwd=2)
arrows(DUEM.reg, F_d.K21b.reg+F_d.K21b.reg.negsd, DUEM.reg, F_d.K21b.reg+F_d.K21b.reg.possd, length=0.1, angle=90, code=3, col='darkorange4', lwd=2)
#arrows(F_d.K14.reg, F_d.K21b.reg+F_d.K21b.reg.negsd, F_d.K14.reg, F_d.K21b.reg+F_d.K21b.reg.possd, length=0.1, angle=90, code=3, col='red', lwd=2)

print(name.region.Ridley.source.K21)
cor(DUEM.reg, F_d.K21b.reg)^2; rmse(DUEM.reg, F_d.K21b.reg)
#cor(F_d.K14.reg, F_d.K21b.reg)^2; rmse(F_d.K14.reg, F_d.K21b.reg)

# correlation between gridded DustCOMM and gridded MERRA-2
cor(mat_2dto1d(DUEM.re), mat_2dto1d(F_d.K21[indi.K,indj.K]), 'complete.obs'); rmse(mat_2dto1d(DUEM.re), mat_2dto1d(F_d.K21[indi.K,indj.K]))
plot.field.log(DUEM.re/1567*5000, lon.K21.r,lat.K21.r, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
plot.field.log(F_d.K21[indi.K,indj.K]/5720*5000, lon.K21.r,lat.K21.r, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
########################################################
# Pu and Ginoux et al. (2020) dust emission threshold

nc = nc_open('/Volumes/GoogleDrive/My Drive/Pu_dust_emis_threshold/Global_wind_erosion_threshold_monthly.nc')
#nc = nc_open('/Volumes/GoogleDrive/My Drive/Pu_dust_emis_threshold/Global_wind_erosion_threshold_annual.nc')

lon_Pu = ncvar_get(nc, 'lon'); lat_Pu = ncvar_get(nc, 'lat')
vv = ncvar_get(nc, 'vv')  # threshold wind speed in m / s
dim(vv)

vv_ann = apply(vv, c(1,2), mean, na.rm=T)
ind_Pu.i = c(22:720,1:21); ind_Pu.j = 65:350
lon_Pu_r = lon_Pu[ind_Pu.i]; lon_Pu_r[700:720] = lon_Pu_r[700:720] + 360
lat_Pu_r = lat_Pu[ind_Pu.j]

plot.field(vv_ann[ind_Pu.i,ind_Pu.j], lon_Pu_r, lat_Pu_r, col=turbo(64))
plot.field(vv_ann[ind_Pu.i,ind_Pu.j], lon_Pu_r, lat_Pu_r)

nc = nc_open('/Volumes/GoogleDrive/My Drive/Pu_dust_emis_threshold/Global_wind_erosion_threshold_monthly_DODtresh05.nc')
#nc = nc_open('/Volumes/GoogleDrive/My Drive/Pu_dust_emis_threshold/Global_wind_erosion_threshold_annual_DODtresh05.nc')

plot.field(vv_ann_05[ind_Pu.i,ind_Pu.j], lon_Pu_r, lat_Pu_r, col=turbo(64))
plot.field(vv_ann_05[ind_Pu.i,ind_Pu.j], lon_Pu_r, lat_Pu_r)
#plot.field(vv[,,3], lon_Pu, lat_Pu, col=turbo(64))
#plot.field(vv[,,9], lon_Pu, lat_Pu, col=turbo(64))

plot.field(vv_ann_05[ind_Pu.i,ind_Pu.j] - vv_ann[ind_Pu.i,ind_Pu.j], lon_Pu_r, lat_Pu_r, col=TEMP_DIFF_65, type='sign', zlim=c(-10,10))

# our scheme's thresholds
# from MERRA2_dustmod_05_06_analysis_new.R 12 Jan 2023
directory = 'L21_final'
filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006_new/', directory, '/thresholds_05x0625_2006_annual_mean.RData', sep='')
load(filename) 
plot.field(ustar_it[i.M2,], lon.M2.r, lat.M2)
plot.field(ustar_ft.wet[i.M2,], lon.M2.r, lat.M2)
plot.field((ustar_ft.wet/ustar_it)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(1,4.5))
plot.field((ustar_ft.wet/ustar_it)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(1,4.5), col=CBR_DRYWET_65[c(4:61)])
# move to saltation level of 0.1 m
plot.field(ustar_it[i.M2,] / 0.4 * log(0.1/1e-4) , lon.M2.r, lat.M2)
plot.field(ustar_ft.wet[i.M2,] / 0.4 * log(0.1/1e-4) , lon.M2.r, lat.M2)

# get annual mean hybrid F_eff
# from MERRA2_dust_05_06_exec.R 12 Jan 2023
filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006_new/', directory, '/F_eff_hybrid_05x0625_2006_annual_mean.RData', sep='')
#save('SSR', 'SSR.1', 'LAI', 'LAI.1', 'F_eff.LC', 'F_eff.LC.1' , 'lon.M2', 'lat.M2', file=filename)
load(filename)
# plot impact threshold modified by hybrid drag partition effect
# move to 10-m height
plot.field(ustar_it[i.M2,] / 0.4 * log(0.1/1e-4), lon.M2.r, lat.M2)
F_eff.LC.1.no0 = F_eff.LC.1; F_eff.LC.1.no0[which(F_eff.LC.1.no0 < 0.3)] = NaN
plot.field((ustar_it / 0.4 * log(10/1e-4) / F_eff.LC.1.no0)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(1,15))

# regrid Pu's threshold to make it comparable with ours
vv_re = sp.dissolve(vv_ann[ind_Pu.i,ind_Pu.j], lon_Pu_r, lat_Pu_r, lon.M2.r, lat.M2)
vv_05_re = sp.dissolve(vv_ann_05[ind_Pu.i,ind_Pu.j], lon_Pu_r, lat_Pu_r, lon.M2.r, lat.M2)
#plot.field(vv_re, lon.M2.r, lat.M2, type='def', zlim=c(1,15))

# ratio between Pu's threshold and our threshold
plot.field(vv_re, lon.M2.r, lat.M2, type='def', zlim=c(1,15))
plot.field(vv_re / (ustar_it / 0.4 * log(10/1e-4) / F_eff.LC.1.no0)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(0,2), col=TEMP_DIFF_65) # in ratio
plot.field.log(vv_re / (ustar_it / 0.4 * log(10/1e-4) / F_eff.LC.1.no0)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-0.478,0.478), col=TEMP_DIFF_65,ticks=c(0.33,0.5,1,2,3), mgp=c(0.4,0.0,0), legend.mar=4) # in ratio (log scale)
plot.field.log(vv_05_re / (ustar_it / 0.4 * log(10/1e-4) / F_eff.LC.1.no0)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-0.478,0.478), col=TEMP_DIFF_65,ticks=c(0.33,0.5,1,2,3), mgp=c(0.4,0.0,0), legend.mar=4) # in ratio (log scale)
plot.field(100*(vv_re / (ustar_it / 0.4 * log(10/1e-4) / F_eff.LC.1.no0)[i.M2,] - 1), lon.M2.r, lat.M2, type='def', zlim=c(-100,100), col=TEMP_DIFF_65) # in percentage change

# difference between Pu's threshold and our threshold
plot.field(vv_re - (ustar_it / 0.4 * log(10/1e-4) / F_eff.LC.1.no0)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-10,10), col=TEMP_DIFF_65)
plot.field(vv_05_re - (ustar_it / 0.4 * log(10/1e-4) / F_eff.LC.1.no0)[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-10,10), col=TEMP_DIFF_65)
#plot.field((vv_re - (ustar_it / 0.4 * log(0.1/1e-4) / F_eff.LC.1.no0)[i.M2,]) * 0.4 / log(0.1/1e-4), lon.M2.r, lat.M2, type='def', zlim=c(-0.5,0.5), col=TEMP_DIFF_65) # translating back to velocity scale instead of velocity

# difference between Pu's threshold and our fluid threshold
plot.field(vv_re - (ustar_ft.wet / 0.4 * log(0.1/1e-4) )[i.M2,], lon.M2.r, lat.M2, type='def', zlim=c(-10,10), col=TEMP_DIFF_65)
#plot.field((vv_re - (ustar_it / 0.4 * log(0.1/1e-4) / F_eff.LC.1.no0)[i.M2,]) * 0.4 / log(0.1/1e-4), lon.M2.r, lat.M2, type='def', zlim=c(-0.5,0.5), col=TEMP_DIFF_65) # translating back to velocity scale instead of velocity
# scatterplot of Pu's threshold versus our threshold
plot(vv_re, (ustar_it / 0.4 * log(0.1/1e-4) / F_eff.LC.1.no0)[i.M2,], 'p')


################################################################################
# relative importance of each added term in explaining DustCOMM as a benchmark

.env = new.env()
load('MERRA2_dustmod_emissions_05x06_c16Jan22.RData',.env); ls(.env)
F_d.K14_def = .env$F_d.K14_def
F_d.K14_130 = .env$F_d.K14_130
F_d.K14_uft_F_eff = .env$F_d.K14_uft_F_eff
F_d.K14_F_eff = .env$F_d.K14_F_eff
F_d.L21 = .env$F_d.L21

# annual mean LAI (use get.merra2.field from MERRA2_dustmod_05x06_analysis_new.R)
FF = get.merra2.field(files_lnd, date_vec, dt, var='LAI')
LAI.ann = FF$VAR; dim(LAI.ann)
LAI.filt.arid = LAI.ann; LAI.filt.arid[which(LAI.filt.arid>1)] = NaN; LAI.filt.arid[which(LAI.filt.arid <= 1)] = 1
LAI.filt.nonarid = LAI.ann; LAI.filt.nonarid[which(LAI.filt.nonarid <= 1)] = NaN; LAI.filt.nonarid[which(LAI.filt.nonarid > 1)] = 1

# K14
sum(F_d.K14_def*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
sum(F_d.K14_def*LAI.filt.nonarid*area.M2, na.rm=TRUE)/1e9   # nonarid Tg/yr-1
sum(F_d.K14_def*LAI.filt.arid*area.M2, na.rm=TRUE)/1e9   # arid Tg/yr-1
plot.field.log(F_d.K14_def[i.M2,]/29254*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR, ps=15)
F_d.K14_def_re = sp.dissolve(F_d.K14_def[i.M2,], lon.M2.r, lat.M2, lon.K21.r, lat.K21.r)   # regrid
sum(F_d.K14_def_re*area.K21.r, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.K14_def_re/29335*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

# K14 plus Dp = 127 um
sum(F_d.K14_130*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.K14_130[i.M2,]/23984*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR, ps=15)
F_d.K14_130_re = sp.dissolve(F_d.K14_130[i.M2,], lon.M2.r, lat.M2, lon.K21.r, lat.K21.r)   # regrid
sum(F_d.K14_130_re*area.K21.r, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.K14_130_re/24050*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

# K14 plus Dp = 127 um + F_eff
sum(F_d.K14_uft_F_eff*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.K14_uft_F_eff[i.M2,]/2886*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR, ps=15)
F_d.K14_uft_F_eff_re = sp.dissolve(F_d.K14_uft_F_eff[i.M2,], lon.M2.r, lat.M2, lon.K21.r, lat.K21.r)   # regrid
sum(F_d.K14_uft_F_eff_re*area.K21.r, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.K14_uft_F_eff_re/2895*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

# K14 plus Dp = 127 um + F_eff + u_it
sum(F_d.K14_F_eff*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.K14_F_eff[i.M2,]/13181*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR, ps=15)
F_d.K14_F_eff_re = sp.dissolve(F_d.K14_F_eff[i.M2,], lon.M2.r, lat.M2, lon.K21.r, lat.K21.r)   # regrid
sum(F_d.K14_F_eff_re*area.K21.r, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.K14_F_eff_re/13219*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

# K14 plus Dp = 127 um + F_eff + u_it + eta
sum(F_d.L21*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
sum(F_d.L21*LAI.filt.arid*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
sum(F_d.L21*LAI.filt.nonarid*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.L21[i.M2,]/11748*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR, ps=15)
F_d.L21_re = sp.dissolve(F_d.L21[i.M2,], lon.M2.r, lat.M2, lon.K21.r, lat.K21.r)   # regrid
sum(F_d.L21_re*area.K21.r, na.rm=TRUE)/1e9   # Tg/yr-1
plot.field.log(F_d.L21_re/11782*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)

##########################################################
# each experiment regressing DustCOMM

# K21 emissions (copied from above)
plot.field.log(F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-3,0), legend.mar=5, col=turbo(64))

#K14 - K21
plot.field(F_d.K14_def_re/29335*5000 - F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-1,1), legend.mar=3, col=TEMP_DIFF_65)
#plot.field.log.diff(F_d.K14_def_re/29335*5000, F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, legend.mar=4.5, col=TEMP_DIFF_65[9:57], oom=1e-3)
plot.field.log.diff(F_d.K14_def_re/29335*5000, F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, legend.mar=4.5, col=TEMP_DIFF_65[9:57], oom=1e-4)
mean((F_d.K14_def_re/29335*5000 - F_d.K21[indi.K,indj.K]/5693*5000)^2)  # RMSE

#K14+Dp=130 - K21
plot.field(F_d.K14_130_re/24050*5000 - F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-1,1), legend.mar=3, col=TEMP_DIFF_65)
plot.field.log.diff(F_d.K14_130_re/24050*5000, F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, legend.mar=4.5, col=TEMP_DIFF_65[9:57], oom=1e-4)
mean((F_d.K14_130_re/24050*5000 - F_d.K21[indi.K,indj.K]/5693*5000)^2) # RMSE

#K14+Dp=130+F_eff - K21
plot.field(F_d.K14_uft_F_eff_re/2895*5000 - F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-1,1), legend.mar=3, col=TEMP_DIFF_65)
plot.field.log.diff(F_d.K14_uft_F_eff_re/2895*5000, F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, legend.mar=4.5, col=TEMP_DIFF_65[9:57], oom=1e-4)
mean((F_d.K14_uft_F_eff_re/2895*5000 - F_d.K21[indi.K,indj.K]/5693*5000)^2)  #RMSE

#K14+Dp=130+F_eff+u_it - K21
plot.field(F_d.K14_F_eff_re/13219*5000 - F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-1,1), legend.mar=3, col=TEMP_DIFF_65)
plot.field.log.diff(F_d.K14_F_eff_re/13219*5000, F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, legend.mar=4.5, col=TEMP_DIFF_65[9:57], oom=1e-4)
mean((F_d.K14_F_eff_re/13219*5000 - F_d.K21[indi.K,indj.K]/5693*5000)^2) #RMSE

#K14+Dp=130+F_eff+u_it+eta - K21
plot.field(F_d.L21_re/11782*5000 - F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, type='def', zlim=c(-1,1), legend.mar=3, col=TEMP_DIFF_65)
plot.field.log.diff(F_d.L21_re/11782*5000, F_d.K21[indi.K,indj.K]/5693*5000, lon.K21.r, lat.K21.r, legend.mar=4.5, col=TEMP_DIFF_65[9:57], oom=1e-4)
mean((F_d.L21_re/11782*5000 - F_d.K21[indi.K,indj.K]/5693*5000)^2)  #RMSE


# spatial variability of each added effect in DustCOMM resolution
eff_Dp130 = F_d.K14_130_re/24050*5000 - F_d.K14_def_re/29335*5000; range(eff_Dp130)
plot.field.log.diff(F_d.K14_130_re/24050*5000, F_d.K14_def_re/29335*5000, lon.K21.r, lat.K21.r, legend.mar=4.5, col=TEMP_DIFF_65[9:57], oom=1e-4)
eff_F_eff = F_d.K14_uft_F_eff_re/2895*5000 - F_d.K14_130_re/24050*5000; range(eff_F_eff)
plot.field.log.diff(F_d.K14_uft_F_eff_re/2895*5000, F_d.K14_130_re/24050*5000, lon.K21.r, lat.K21.r, legend.mar=4.5, col=TEMP_DIFF_65[9:57], oom=1e-4)
eff_uit = F_d.K14_F_eff_re/13219*5000 - F_d.K14_uft_F_eff_re/2895*5000
plot.field.log.diff(F_d.K14_F_eff_re/13219*5000, F_d.K14_uft_F_eff_re/2895*5000, lon.K21.r, lat.K21.r, legend.mar=4.5, col=TEMP_DIFF_65[9:57], oom=1e-4)
eff_int = F_d.L21_re/11782*5000 - F_d.K14_F_eff_re/13219*5000; range(eff_int)
plot.field.log.diff(F_d.L21_re/11782*5000, F_d.K14_F_eff_re/13219*5000, lon.K21.r, lat.K21.r, legend.mar=4.5, col=TEMP_DIFF_65[9:57], oom=1e-4)
eff_int_tot = F_d.L21_re/11782*5000 - F_d.K14_uft_F_eff_re/2895*5000; range(eff_int_tot)
plot.field.log.diff(F_d.L21_re/11782*5000, F_d.K14_uft_F_eff_re/2895*5000, lon.K21.r, lat.K21.r, legend.mar=4.5, col=TEMP_DIFF_65[9:57], oom=1e-4)

################################################################################
# spatial variability of each added effect in MERRA-2 resolution (normalized)

# annual mean LAI (use get.merra2.field from MERRA2_dustmod_05x06_analysis_new.R)
FF = get.merra2.field(files_lnd, date_vec, dt, var='LAI')
LAI.ann = FF$VAR; dim(LAI.ann)
LAI.filt.arid = LAI.ann; LAI.filt.arid[which(LAI.filt.arid>1)] = NaN; LAI.filt.arid[which(LAI.filt.arid <= 1)] = 1
LAI.filt.nonarid = LAI.ann; LAI.filt.nonarid[which(LAI.filt.nonarid <= 1)] = NaN; LAI.filt.nonarid[which(LAI.filt.nonarid > 1)] = 1

eff_Dp130 = F_d.K14_130/24050*5000 - F_d.K14_def/29335*5000; range(eff_Dp130)
eff_Dp130_no0 = eff_Dp130; eff_Dp130_no0[which(eff_Dp130_no0==0)] = NaN; sd(eff_Dp130_no0*area.M2/1e9, na.rm=T)
plot.field.log.diff((F_d.K14_130/24050*5000)[i.M2,], (F_d.K14_def/29335*5000)[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=0.01); 
sum(abs(eff_Dp130)*area.M2/1e9) # global
sum(abs(eff_Dp130*LAI.filt.arid)*area.M2/1e9, na.rm=T) # arid
sum(abs(eff_Dp130*LAI.filt.nonarid)*area.M2/1e9, na.rm=T) # nonarid

eff_F_eff = F_d.K14_uft_F_eff/2895*5000 - F_d.K14_130/24050*5000; range(eff_F_eff)
eff_F_eff_no0 = eff_F_eff; eff_F_eff_no0[which(eff_F_eff_no0==0)] = NaN; sd(eff_F_eff_no0*area.M2/1e9, na.rm=T)
plot.field.log.diff(F_d.K14_uft_F_eff[i.M2,]/2895*5000, F_d.K14_130[i.M2,]/24050*5000, lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=1); 
sum(abs(eff_F_eff)*area.M2/1e9) # global
sum(abs(eff_F_eff*LAI.filt.arid)*area.M2/1e9, na.rm=T) # arid
sum(abs(eff_F_eff*LAI.filt.nonarid)*area.M2/1e9, na.rm=T) # nonarid

eff_uit = F_d.K14_F_eff/13219*5000 - F_d.K14_uft_F_eff/2895*5000
eff_uit_no0 = eff_uit; eff_uit_no0[which(eff_uit==0)] = NaN; sd(eff_uit*area.M2/1e9, na.rm=T)
plot.field.log.diff(F_d.K14_F_eff[i.M2,]/13219*5000, F_d.K14_uft_F_eff[i.M2,]/2895*5000, lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=1); sum(abs(eff_uit)*area.M2/1e9)

eff_int = F_d.L21/11782*5000 - F_d.K14_F_eff/13219*5000; range(eff_int)
eff_int_no0 = eff_int; eff_int_no0[which(eff_int==0)] = NaN; sd(eff_int*area.M2/1e9, na.rm=T)
plot.field.log.diff(F_d.L21[i.M2,]/11782*5000, F_d.K14_F_eff[i.M2,]/13219*5000, lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=0.01); sum(abs(eff_int)*area.M2/1e9)

eff_int_tot = F_d.L21/11782*5000 - F_d.K14_uft_F_eff/2895*5000; range(eff_int_tot)
eff_int_tot_no0 = eff_int_tot; eff_int_tot_no0[which(eff_int_tot==0)] = NaN; sd(eff_int_tot*area.M2/1e9, na.rm=T)
plot.field.log.diff(F_d.L21[i.M2,]/11782*5000, F_d.K14_uft_F_eff[i.M2,]/2895*5000, lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=1); 
sum(abs(eff_int_tot)*area.M2/1e9) # global
sum(abs(eff_int_tot*LAI.filt.arid)*area.M2/1e9, na.rm=T) # arid
sum(abs(eff_int_tot*LAI.filt.nonarid)*area.M2/1e9, na.rm=T) # nonarid

eff_L21 = F_d.L21/11782*5000 - F_d.K14_def/29335*5000; range(eff_int_tot)
eff_L21_no0 = eff_L21; eff_L21_no0[which(eff_L21==0)] = NaN; sd(eff_L21*area.M2/1e9, na.rm=T)
plot.field.log.diff(F_d.L21[i.M2,]/11782*5000, F_d.K14_def[i.M2,]/29335*5000, lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=1); 
sum(abs(eff_L21)*area.M2/1e9) # global
sum(abs(eff_L21*LAI.filt.arid)*area.M2/1e9, na.rm=T) # arid
sum(abs(eff_L21*LAI.filt.nonarid)*area.M2/1e9, na.rm=T) # nonarid

##################################################################################
# spatial variability of each added effect in MERRA-2 resolution (not normalized)

eff_Dp130 = F_d.K14_130 - F_d.K14_def; range(eff_Dp130)
plot.field.log.diff((F_d.K14_130)[i.M2,], (F_d.K14_def)[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=10)
eff_F_eff = F_d.K14_uft_F_eff - F_d.K14_130; range(eff_F_eff)
plot.field.log.diff(F_d.K14_uft_F_eff[i.M2,], F_d.K14_130[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=10)
eff_uit = F_d.K14_F_eff - F_d.K14_uft_F_eff
plot.field.log.diff(F_d.K14_F_eff[i.M2,], F_d.K14_uft_F_eff[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=10)
eff_int = F_d.L21 - F_d.K14_F_eff; range(eff_int)
plot.field.log.diff(F_d.L21[i.M2,], F_d.K14_F_eff[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=10)
eff_int_tot = F_d.L21 - F_d.K14_uft_F_eff; range(eff_int_tot)
plot.field.log.diff(F_d.L21[i.M2,], F_d.K14_uft_F_eff[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=10)

eff_L21 = F_d.L21 - F_d.K14_def; range(eff_int_tot)
plot.field.log.diff(F_d.L21[i.M2,], F_d.K14_def[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], oom=1e-5, oom.max=10)
plot.field(F_d.L21[i.M2,] - F_d.K14_def[i.M2,], lon.M2.r, lat.M2, legend.mar=5, col=TEMP_DIFF_65[9:57], type='sign', zlim=c(-2,2))

###################################################################
# in making Table for arid and nonarid emissions, refer to 
# MERRA2_dustmod_05x06_analysis.R to use get.tot.dustemis with DPM='Z03'



###################################################################
# making Table for emissions from each region

# function for getting annual total dust emissions
get.tot.dustemis = function(date_vec, dt, DPM=DPM, use.intermittency = NULL) {
  # get annual dust emission total per grid (kg m-2 yr-1) 
  .env1 = new.env()
  for (d in 1:length(date_vec)) {
    print(date_vec[d])
    if (DPM=='K14') {filename = paste('/Volumes/SEAGATE/MERRA2/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_globe_fragexpt=',frag_expt_lim,'_intmtncy=',use.intermittency,'/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')}
    else if (DPM=='Z03') {filename = paste('/Volumes/GoogleDrive/My Drive/upscale_2006/OUTPUTS_',clay_dataset,'_',LC,'_',DPM,'_src_fn=',src_fn,'/dust_emis_flx_05x0625_', date_vec[d], '.RData', sep='')}
    load(filename, envir=.env1)
    #if (d==1) {F_d = .env1$F_d}
    #else if (d >= 2) {F_d = F_d + .env1$F_d}
    if (d==1) {F_d=.env1$F_d; F_d[which(is.na(F_d))]=0}
    else if (d >= 2) {F_d1=.env1$F_d; F_d1[which(is.na(F_d1))]=0; F_d=F_d+F_d1}
    
    #F_d = abind(F_d, (apply(.env1$F_d, c(1,2), sum, na.rm=T)*3600*dt), along=3)
    #F_d = apply(F_d, c(1,2), sum, na.rm=T)
  }
  dim(F_d)
  F_d = apply(F_d, c(1,2), sum, na.rm=T)*3600*dt
  print(paste('total global dust emission is ',sum(F_d*area.M2, na.rm=TRUE)/1e9, ' Tg/yr',sep=''))
  return(list(F_d=F_d,lon.M2=lon.M2, lat.M2=lat.M2))
}


# plot emissions for Z03-G scheme 
src_fn='Ginoux'; clay_dataset = 'FAO'
FF = get.tot.dustemis(date_vec, dt, DPM='Z03')  # this is Z03-G scheme
F_d.Z03.G = FF$F_d
sum(F_d.Z03.G*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
sum(F_d.Z03.G*LAI.filt.arid*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1, arid region
sum(F_d.Z03.G*LAI.filt.nonarid*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1, nonarid region
plot.field.log(F_d.Z03.G[i.M2,]/442*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR, ps=15)
# Patagonia emissions
Z03.G.Pat = extract.value.selfdef.region(F_d.Z03.G[i.M2,]/442*5000, lon.M2.r, lat.M2, low.lat=(-58), high.lat=(-38), left.lon=(-80), right.lon=(-65), output.matrix=TRUE)
Z03.G.Pat = sum(Z03.G.Pat$matrix.2d * area.M2/1e9, na.rm=T)  # total Patagonia emissions in Tg / yr
Z03.G.region = extract.value.region(F_d.Z03.G[i.M2,]/442*5000, lon.M2.r, lat.M2, output.matrix.3d=TRUE, frame='source')
names(Z03.G.region)
Z03.G.region.grid = Z03.G.region$matrix.3d  # gridded emissions separated by regions
Z03.G.reg = apply(sweep(Z03.G.region.grid, c(1,2), area.M2/1e9, '*'), c(3), sum, na.rm=T)
Z03.G.reg[10] = Z03.G.reg[10] + Z03.G.Pat
# percentages for Z03-G emissions
Z03.G.reg / sum(Z03.G.reg[1:9]) * 100


# plot emissions for Z03-G scheme 
F_d.Z03.Z = F_d.Z03.G/mbl_bsn_fct.g*mbl_bsn_fct.z
sum(F_d.Z03.Z*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1
sum(F_d.Z03.Z*LAI.filt.arid*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1, arid
sum(F_d.Z03.Z*LAI.filt.nonarid*area.M2, na.rm=TRUE)/1e9   # Tg/yr-1, nonarid
plot.field.log(F_d.Z03.Z[i.M2,]/423*5000, lon.M2.r, lat.M2, type='def', zlim=c(-4,0), legend.mar=5, col=WBGYR)
# Patagonia emissions
Z03.Z.Pat = extract.value.selfdef.region(F_d.Z03.Z[i.M2,]/442*5000, lon.M2.r, lat.M2, low.lat=(-58), high.lat=(-38), left.lon=(-80), right.lon=(-65), output.matrix=TRUE)
Z03.Z.Pat = sum(Z03.Z.Pat$matrix.2d * area.M2/1e9, na.rm=T)  # total Patagonia emissions in Tg / yr
Z03.Z.region = extract.value.region(F_d.Z03.Z[i.M2,]/442*5000, lon.M2.r, lat.M2, output.matrix.3d=TRUE, frame='source')
names(Z03.Z.region)
Z03.Z.region.grid = Z03.Z.region$matrix.3d  # gridded emissions separated by regions
Z03.Z.reg = apply(sweep(Z03.Z.region.grid, c(1,2), area.M2/1e9, '*'), c(3), sum, na.rm=T)
Z03.Z.reg[10] = Z03.Z.reg[10] + Z03.Z.Pat
# percentages for Z03-Z emissions
Z03.Z.reg / sum(Z03.Z.reg[1:9]) * 100

# get regional emissions for our scheme 
# Patagonia emissions
F_d.L21.Pat = extract.value.selfdef.region(F_d.L21[i.M2,]/11748*5000, lon.M2.r, lat.M2, low.lat=(-58), high.lat=(-38), left.lon=(-80), right.lon=(-65), output.matrix=TRUE)
F_d.L21.Pat = sum(F_d.L21.Pat$matrix.2d * area.M2/1e9, na.rm=T)  # total Patagonia emissions in Tg / yr
#plot.field.log(F_d.L21.Pat$matrix.2d, lon.M2.r, lat.M2, legend.mar=5)
L21.region = extract.value.region(F_d.L21[i.M2,]/11748*5000, lon.M2.r, lat.M2, output.matrix.3d=TRUE, frame='source')
names(L21.region)
L21.region.grid = L21.region$matrix.3d  # gridded emissions separated by regions
F_d.L21.reg = apply(sweep(L21.region.grid, c(1,2), area.M2/1e9, '*'), c(3), sum, na.rm=T)
F_d.L21.reg[10] = F_d.L21.reg[10] + F_d.L21.Pat
# percentages for regional emissions
F_d.L21.reg / sum(F_d.L21.reg[1:9]) * 100


K14.region = extract.value.region(F_d.K14_def[i.M2,]/29254*5000, lon.M2.r, lat.M2, output.matrix.3d=TRUE, frame='source')
K14.region.grid = K14.region$matrix.3d  # gridded emissions separated by regions
plot.field.log(K14.region.grid[,,1], lon.M2.r, lat.M2, type='def', zlim=c(-3,0), legend.mar=5, col=WBGYR)
F_d.K14.reg = apply(sweep(K14.region.grid, c(1,2), area.M2/1e9, '*'), c(3), sum, na.rm=T)
# Patagonia emissions
F_d.K14.Pat = extract.value.selfdef.region(F_d.K14_def[i.M2,]/29254*5000, lon.M2.r, lat.M2, low.lat=(-58), high.lat=(-39), left.lon=(-80), right.lon=(-65), output.matrix=TRUE)
F_d.K14.Pat = sum(F_d.K14.Pat$matrix.2d * area.M2/1e9, na.rm=T)  # total Patagonia emissions in Tg / yr
F_d.K14.reg[10] = F_d.K14.reg[10] + F_d.K14.Pat
# percentages for regional emissions
F_d.K14.reg / sum(F_d.K14.reg[1:9]) * 100




# percentages for MERRA-2 dust emissions
DUEM.reg / sum(DUEM.reg[1:9]) * 100
