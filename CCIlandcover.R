
setwd('/Volumes/GoogleDrive/My Drive/CCI_Landcover')
#library(doParallel); library(foreach); library(parallel)
library(raster)
library(ncdf4)
library(fields); library(maps)
source('/Volumes/GoogleDrive/My Drive/CESM/CESMDustinR/get_geo.R')
source('/Volumes/GoogleDrive/My Drive/CESM/CESMDustinR/sptial_plot_fns.R')
load('/Volumes/GoogleDrive/My Drive/CESM/CESMDustinR/WBGYR_scheme.RData'); load('/Volumes/GoogleDrive/My Drive/CESM/CESMDustinR/TEMP_DIFF.RData')


nc = nc_open('ESACCI-LC-L4-LCCS-Map-300m-P1Y-2006-v2.0.7b.nc')

print(nc)
names(nc$var$lccs_class)



class_values = c(0,10,11,12,20,30,40,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,-126,-116,-106,-105,-104,-103,-96,-86,-76,-66,-56,-55,-54,-46,-36)
length(class_values)

class_meanings = c('no_data', 'cropland_rainfed', 'cropland_rainfed_herbaceous_cover', 'cropland_rainfed_tree_or_shrub_cover', 'cropland_irrigated', 'mosaic_cropland', 'mosaic_natural_vegetation', 'tree_broadleaved_evergreen_closed_to_open', 'tree_broadleaved_deciduous_closed_to_open', 'tree_broadleaved_deciduous_closed', 'tree_broadleaved_deciduous_open', 'tree_needleleaved_evergreen_closed_to_open', 'tree_needleleaved_evergreen_closed', 'tree_needleleaved_evergreen_open', 'tree_needleleaved_deciduous_closed_to_open', 'tree_needleleaved_deciduous_closed', 'tree_needleleaved_deciduous_open', 'tree_mixed', 'mosaic_tree_and_shrub', 'mosaic_herbaceous', 'shrubland', 'shrubland_evergreen', 'shrubland_deciduous', 'grassland', 'lichens_and_mosses', 'sparse_vegetation', 'sparse_tree', 'sparse_shrub', 'sparse_herbaceous', 'tree_cover_flooded_fresh_or_brakish_water', 'tree_cover_flooded_saline_water', 'shrub_or_herbaceous_cover_flooded', 'urban', 'bare_areas', 'bare_areas_consolidated', 'bare_areas_unconsolidated', 'water', 'snow_and_ice')



length(ncvar_get(nc, 'lon'))/40
length(ncvar_get(nc, 'lat'))/40

start = c(1,1)
count = c(3240,1620)

ptm = proc.time()
lccs_class = ncvar_get(nc, 'lccs_class', start=start, count=count)
proc.time() - ptm

lon.cci = ncvar_get(nc, 'lon')#[start[1]+c(1:3240)-1]
lat.cci = ncvar_get(nc, 'lat')#[start[1]+c(1:1620)-1]

plot.field(lccs_class, lon.cci, lat.cci)


indi = 1:576; indj = 65:350 
# get lon lat
nc1 = nc_open('/Volumes/GoogleDrive/My Drive//MERRA2/upscale_2006/MERRA2_300.tavg1_2d_lnd_Nx.20060101.SUB.nc')
lon.M2 = ncvar_get(nc1, 'lon')[indi]; lat.M2 = ncvar_get(nc1, 'lat')[indj]
nc_close(nc1)
i.M2 = c(20:576,1:19)
lon.M2.r = lon.M2[i.M2]; lon.M2.r[558:576]=lon.M2.r[558:576]+360


dlon = lon.M2[2]-lon.M2[1]
dlat = lat.M2[2]-lat.M2[1]

#bare_u = bare_c = bare = water = array(NaN, dim=c(length(lon.M2), length(lat.M2)))
dim(bare)

landcover = array(NaN, dim=c(length(class_values), length(lon.M2), length(lat.M2)))


ptm = proc.time()
for (i in 1:length(lon.M2)) {
  print(lon.M2[i])
  for (j in 1:length(lat.M2)) {
#for (i in 1:10) {
#  for (j in 1:10) {    
    indi = which(lon.cci >= lon.M2[i]-dlon/2 & lon.cci <= lon.M2[i]+dlon/2)
    indj = which(lat.cci >= lat.M2[j]-dlat/2 & lat.cci <= lat.M2[j]+dlat/2)
    
    start = c(min(indi), min(indj))
    count = c(max(indi)-min(indi)+1, max(indj)-min(indj)+1)
    #dim(lccs_class)
    lccs_class = ncvar_get(nc, 'lccs_class', start=start, count=count)
    #water[i,j] = length(which(lccs_class == (-46))) / length(lccs_class)
    
    tot.grid = length(lccs_class)
    for (k in 1:length(class_values)) {
      
      landcover[k,i,j] = length(which(lccs_class == class_values[k])) / tot.grid
      
    }
    
    #bare[i,j] = length(which(lccs_class == (-56))) / length(lccs_class)
    #bare_c[i,j] = length(which(lccs_class == (-55))) / length(lccs_class)
    #bare_u[i,j] = length(which(lccs_class == (-54))) / length(lccs_class)
  }
}
proc.time() - ptm

#plot.field(bare, lon.M2, lat.M2)
#plot.field(bare_c, lon.M2, lat.M2)
#plot.field(water, lon.M2, lat.M2)

# bare land
plot.field(landcover[34,,], lon.M2, lat.M2, col=WBGYR) # bare
plot.field(landcover[35,,], lon.M2, lat.M2, col=WBGYR) # bare consolidated
plot.field(landcover[36,,], lon.M2, lat.M2, col=WBGYR) # bare unconsolidated
land_bare = landcover[34,,]+landcover[35,,]+landcover[36,,]
plot.field(land_bare[i.M2,], lon.M2.r, lat.M2, col=WBGYR) # bare total


# short plants
class_meanings

# sparse plants
land_sparse = landcover[26,,]+landcover[27,,]+landcover[28,,]+landcover[29,,]
plot.field(land_sparse[i.M2,], lon.M2.r, lat.M2, col=WBGYR) 

# cropland
land_crop = landcover[2,,]+landcover[3,,]+landcover[4,,]+landcover[5,,]
plot.field(land_crop[i.M2,], lon.M2.r, lat.M2, col=WBGYR) 

# mosaic
land_mosaic = landcover[6,,]+landcover[7,,]+landcover[19,,]+landcover[20,,]
plot.field(land_mosaic[i.M2,], lon.M2.r, lat.M2, col=WBGYR) 
# lichens/mosses
land_lichensmosses = landcover[25,,]
plot.field(land_lichensmosses[i.M2,], lon.M2.r, lat.M2, col=WBGYR) 

# grassland
land_grassland = landcover[24,,]
plot.field(land_grassland[i.M2,], lon.M2.r, lat.M2, col=WBGYR) 

# shrubland
land_shrub = landcover[21,,]+landcover[22,,]+landcover[23,,]
plot.field(land_shrub[i.M2,], lon.M2.r, lat.M2, col=WBGYR) 

# tree_mixed
land_mixedtree = landcover[18,,]
plot.field(land_mixedtree[i.M2,], lon.M2.r, lat.M2, col=WBGYR)

# shrub_or_herbaceous_cover_flooded
land_floodshruborherbaceous = landcover[32,,]
plot.field(land_floodshruborherbaceous[i.M2,], lon.M2.r, lat.M2, col=WBGYR)

# snow_and_ice
land_snowice = landcover[38,,]
plot.field(land_snowice[i.M2,], lon.M2.r, lat.M2, col=WBGYR)

# tree_cover_flooded
land_floodtree = landcover[30,,]+landcover[31,,]
plot.field(land_floodtree[i.M2,], lon.M2.r, lat.M2, col=WBGYR)

# needleleaved evergreen
land_needlegreen = landcover[12,,]+landcover[13,,]+landcover[14,,]
plot.field(land_needlegreen[i.M2,], lon.M2.r, lat.M2, col=WBGYR)

# needleleaved deciduous
land_needledeciduous = landcover[15,,]+landcover[16,,]+landcover[17,,]
plot.field(land_needledeciduous[i.M2,], lon.M2.r, lat.M2, col=WBGYR)

# needleleaved open
#land_needleopen = landcover[12,,]+landcover[14,,]+landcover[15,,]+landcover[17,,]
land_needleopen = landcover[15,,]
#land_needleopen = landcover[12,,]+landcover[15,,]
plot.field(land_needleopen[i.M2,], lon.M2.r, lat.M2, col=WBGYR)

# all short vegetation
land_shortveg = land_crop + land_grassland + land_mosaic + land_mixedtree + land_floodshruborherbaceous + land_shrub + land_sparse + land_lichensmosses + land_needleopen
plot.field(land_shortveg[i.M2,], lon.M2.r, lat.M2, col=WBGYR)


save('landcover','lon.M2','lat.M2',file='CCI_frc_area.RData')

.env = new.env()
load('CCI_frc_area.RData', envir=.env)

ls(.env)
lon.M2 = .env$lon.M2
lat.M2 = .env$lat.M2
landcover = .env$landcover[34,,]

plot.field(landcover, lon.M2, lat.M2)
plot.field.proj(landcover, lon.M2, lat.M2)
