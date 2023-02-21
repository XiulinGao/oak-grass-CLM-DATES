## CA-OakWoodland-domain.R

## create a mask for oak woodland and remap to WRF domain

library(sp)
library(sf)
library(raster)
library(ncdf4)
library(tidyverse)
library(stars)
library(lubridate)

veg_path = file.path("~/Google Drive/My Drive/CA-OakWoodland/ds2723.gdb")
grass_path = file.path("~/Google Drive/My Drive/CA-OakWoodland/wrf-extent-masked-grassland.tif") #mask created from mask-wrf-domain.R using NLCD data
wrf_path = file.path("~/Google Drive/My Drive/9km-WRF-1980-2020/1981-01.nc")
wgs      = "+init=EPSG:4326"
my_theme = theme_bw() + theme(panel.ontop=TRUE, panel.background=element_blank())

author_name  = "Xiulin Gao"
author_email = "xiulingao@lbl.gov"
undef        = -9999.0000

wrf_t = nc_open(wrf_path)
tobt  = ncvar_get(wrf_t,"TBOT")
XLONG = ncvar_get(wrf_t,"LONGXY")
XLAT  = ncvar_get(wrf_t, "LATIXY")
nc_close(wrf_t)
tobt     = tobt[,,1]
tobt_vec = as.vector(tobt)
x_vec    = as.vector(XLONG)
y_vec    = as.vector(XLAT)
wrf_df   = as.data.frame(cbind(x_vec,y_vec,tobt_vec))
id       = 1:nrow(wrf_df)
wrf_df$cellid = id

layers = st_layers(dsn = veg_path)
oakwl  = st_read(veg_path,layer=layers$name[1])
oakwl  = oakwl %>% filter(OakWoodld=="Y") # only oak woodland
oakwl_wgs = st_transform(oakwl,crs=4326)

wrf_extnt = raster::extent(-130.2749,-108.9862,28.62024,45.71207)
wrf_rs    = raster(wrf_extnt,ncols=147,nrows=151,crs=wgs)
val       = 1:ncell(wrf_rs)
wrf_rs = setValues(wrf_rs,val)
oak_wrf   = mask(wrf_rs,oakwl_wgs)
oak_binary = oak_wrf
oak_binary[oak_binary > 0] = 1
oak_df     = as.data.frame(oak_binary,xy=TRUE)


oak_qry          = wrf_df %>% dplyr::select(x_vec,y_vec,cellid,tobt_vec)
oak_qry$oak_mask = 0

#search for the nearest neighbor in wrf for each corresponding filtered nlcd cell
nn_pt        = RANN::nn2(oak_df[,1:2],oak_qry[,1:2],k=1) 
oak_qry$id   = as.vector(nn_pt$nn.idx) #add resulted row index of wrf to nlcd
oak_qry$dist = as.vector(nn_pt$nn.dists)

oak_qry = oak_qry                                     %>% 
  mutate(if_oak = oak_df$layer[id])                   %>% 
  mutate(oak_mask = ifelse(is.na(if_oak),0,1))                 


wrf_dfil = oak_qry                                     %>%  
  dplyr::select(x_vec,y_vec,oak_mask,tobt_vec)         %>% 
  rename(lon=x_vec,lat=y_vec)                          %>% 
  mutate(tobt_vec = ifelse(oak_mask==0,NA,tobt_vec))

fil_var   = matrix(wrf_dfil$tobt_vec, nrow=147,ncol=151)
x_arr     = matrix(wrf_dfil$lon,nrow=147,ncol=151)
y_arr     = matrix(wrf_dfil$lat,nrow=147,ncol=151)
wrf_star  = st_as_stars(fil_var)
wrf_star  = st_as_stars(wrf_star, curvilinear=list(X1=x_arr,X2=y_arr), crs=wgs)
wrf_sf    = st_as_sf(wrf_star,as_points=FALSE,na.rm=FALSE)
ca_co     = USAboundaries::us_counties(resolution = "high", states = "CA")

##plot to see how the active domain looks like

ggplot() + geom_sf(data=wrf_sf,colour="grey50", aes(fill=A1),lwd=0)+
  coord_sf(crs=st_crs(wgs)) + 
  my_theme +scale_fill_continuous(low="thistle2", high="darkred", 
                                  guide="colorbar",na.value="grey50")+
  geom_sf(data = ca_co, color = alpha("black", alpha=0.2),lwd=0.1,fill=NA) +
  geom_point(aes(x=-120.9508,y=38.4133), colour=alpha("blue",0.6), size=0.2)

## merge oak woodland and annual grassland

grass_domain = raster(grass_path)
grass_df     = as.data.frame(grass_domain,xy=TRUE)

grass_qry        = wrf_df %>% dplyr::select(x_vec,y_vec,cellid,tobt_vec)
grass_qry$mask80 = 0

#search for the nearest neighbor in wrf for each corresponding filtered nlcd cell
nn_pt          = RANN::nn2(grass_df[,1:2],grass_qry[,1:2],k=1) 
grass_qry$id   = as.vector(nn_pt$nn.idx) #add resulted row index of wrf to nlcd
grass_qry$dist = as.vector(nn_pt$nn.dists)

grass_qry = grass_qry                                     %>% 
  mutate(cover = grass_df$layer[id])                      %>% 
  mutate(mask80 = ifelse(is.na(cover),0,1))               %>% 
  mutate(mask80 = ifelse(y_vec>=34 & y_vec<=40 & x_vec>=-124.5 & x_vec<=-114.5
                         ,mask80,0))                      

## merge the two masks

merg_mask = oak_qry %>% dplyr::select(-c(tobt_vec,id,dist)) %>% left_join(grass_qry,by=c("x_vec","y_vec","cellid"))

## as some grid cells are identified as both oak woodland and annual grassland, we 
## re-identify some of the oak cells to be grassland as there are more oak cells

merg_newmask = merg_mask %>% mutate(oak_mask = if_else(mask80>0 , 0, oak_mask),
                                    domain   = if_else(mask80>0 | oak_mask>0, 1, 0))

merg_newmask = merg_newmask %>% dplyr::select(x_vec,y_vec,oak_mask,mask80,domain,cover,if_oak) %>% 
               rename(lon        = x_vec,
                      lat        = y_vec,
                      grass_mask = mask80,
                      grass_cover= cover)

data.table::fwrite(merg_newmask,file = file.path("~/Google Drive/My Drive/CA-OakWoodland/grass-oakwoodland-mask.csv"))

## plot merged domain

plot_df = merg_newmask %>% mutate(type = case_when(grass_mask>0 ~ "Annual grassland",
                                                   oak_mask>0 ~ "Oak woodland",
                                                   grass_mask ==0 & oak_mask==0 ~ "Other"))

wrf_dfil = plot_df %>%  dplyr::select(lon,lat,type)                 

fil_var   = matrix(wrf_dfil$type, nrow=147,ncol=151)
x_arr     = matrix(wrf_dfil$lon,nrow=147,ncol=151)
y_arr     = matrix(wrf_dfil$lat,nrow=147,ncol=151)
wrf_star  = st_as_stars(fil_var)
wrf_star  = st_as_stars(wrf_star, curvilinear=list(X1=x_arr,X2=y_arr), crs=wgs)
wrf_sf    = st_as_sf(wrf_star,as_points=FALSE,na.rm=FALSE)
ca_co     = USAboundaries::us_counties(resolution = "high", states = "CA")

##plot to see how the active domain looks like

ggplot() + geom_sf(data=wrf_sf,colour="grey50", aes(fill=A1),lwd=0)+
  coord_sf(crs=st_crs(wgs)) + 
  my_theme + scale_fill_manual(values = c("Annual grassland" = "blue","Oak woodland"="green","Other"="grey"))+
  geom_sf(data = ca_co, color = alpha("black", alpha=0.2),lwd=0.1,fill=NA) +
  geom_point(aes(x=-120.9508,y=38.4133), colour=alpha("red",0.6), size=0.2)



### the merged domain
land_mask  <- array(merg_newmask$domain,dim=c(147,151))
land_mkdif <- land_mask


## create new nc file as land mask
xx  = ncdim_def( name="lon"   ,units="",vals= sequence(147)  ,create_dimvar=FALSE)
yy  = ncdim_def( name="lat"   ,units="",vals= sequence(151)  ,create_dimvar=FALSE)
nc_xy  = list   (xx,yy)
xy     = c(147,151)
file_name = file.path("~/Google Drive/My Drive/CA-OakWoodland/wrf_CA_grass-oakwoodland.nc")
nc_vlist        = list()
nc_vlist$LONGXY = ncvar_def(  name      = "lsmlon"
                              , units    = "degrees_east"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "longitude"
)#end ncvar_def
nc_vlist$LATIXY = ncvar_def( name       = "lsmlat"
                             , units    = "degrees_north"
                             , dim      = nc_xy
                             , missval  = undef
                             , longname = "latitude"
)#end ncvar_def
nc_vlist$mask1   = ncvar_def( name      = "landmask"
                              , units    = "unitless"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "mask for land domain, 1 being cell active"
)#end ncvar_def
nc_vlist$mask2   = ncvar_def( name      = "mod_lnd_props"
                              , units    = "unitless"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "mask for modifying land property, 1 being active land cell"
)#end ncvar_def

### define global attributes

att_template = list(   title          = "To be replaced when looping through months"
                     , date_created   = paste0(as.character(now(tzone="UTC")), "UTC")
                     , source_code    = "mask-wrf-domain.R"
                     , code_notes     = "land mask for WRF domain created using NLCD herb cover >=80"
                     , code_developer = paste0( author_name
                                                ," <"
                                                , author_email
                                                ,">"
                     )#end paste0
                     , file_author    = paste0(author_name," <",author_email,">")
)#end list

nc_new <- nc_create(filename=file_name,vars=nc_vlist,verbose=FALSE)
dummy = ncvar_put(nc=nc_new,varid="lsmlon",vals=array(data=x_vec,dim=xy))
dummy = ncvar_put(nc=nc_new,varid="lsmlat", vals=array(data=y_vec, dim=xy))
dummy = ncvar_put(nc=nc_new, varid ="landmask",vals=land_mask)
dummy = ncvar_put(nc=nc_new, varid ="mod_lnd_props",vals=land_mkdif)


nc_title   = "Land mask for annual grassland and oakwoodland on WRF domain"
att_global = modifyList( x = att_template, val = list( title = nc_title ))


# Loop through global attributes
for (l in seq_along(att_global)){
  # Current attribute information
  att_name  = names(att_global)[l]
  att_value = att_global[[l]]
  
  # Add attribute 
  dummy = ncatt_put(nc=nc_new,varid=0,attname=att_name,attval=att_value)
}#end for (l in seq_along(att_global))
nc_close(nc_new)
