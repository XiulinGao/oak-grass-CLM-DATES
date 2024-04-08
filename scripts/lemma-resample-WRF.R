## tonzi-benchmark.R
## script for processing benchmark data for tonzi ranch

library(raster)
library(tidyverse)
library(ggplot2)
library(sf)
library(stars)
library(ncdf4)
library(lubridate)

## LEMMA data
obs_main    = file.path("~/Google Drive/My Drive/CA-grassland-simulationDoc/benchmark")
clim_main   = file.path(obs_main,"global-historical-climate/wc2.1_10m_bio")
main_lema   = file.path("~/Google Drive/My Drive/CA-OakWoodland/benchmark/dbf13baf-7da7-49ab-88e0-86dfa968ed2b/rasters/")
wrf_path    = file.path("~/Google Drive/My Drive/9km-WRF-1980-2020/1981-01.nc")
mass_file   = "bphh_ge_3_crm_2017.tif"
cancov_file = "cancov_hdw_2017.tif"
cancovt_file= "cancov_2017.tif"
quagba_file = "quag_ba_2017.tif"
quchba_file = "quch2_ba_2017.tif"
qudoba_file = "qudo_ba_2017.tif"
quenba_file = "quen_ba_2017.tif"
qugaba_file = "quga4_ba_2017.tif"
qukeba_file = "quke_ba_2017.tif"
quloba_file = "qulo_ba_2017.tif"
qumuba_file = "qumu_ba_2017.tif"
quwiba_file = "quwi2_ba_2017.tif"
baall_file  = "ba_ge_3_2017.tif"
bacv_file   = "ba_ge_3_2017_cv.tif"

cancov = raster(file.path(main_lema,cancov_file))
biomas = raster(file.path(main_lema,mass_file))
qugaba = raster(file.path(main_lema,qugaba_file))
quchba = raster(file.path(main_lema,quchba_file))
qudoba = raster(file.path(main_lema,qudoba_file))
quenba = raster(file.path(main_lema,quenba_file))
quagba = raster(file.path(main_lema,quagba_file))
qukeba = raster(file.path(main_lema,qukeba_file))
quloba = raster(file.path(main_lema,quloba_file))
qumuba = raster(file.path(main_lema,qumuba_file))
quwiba = raster(file.path(main_lema,quwiba_file))
ball   = raster(file.path(main_lema,baall_file))
bacv  = raster(file.path(main_lema,bacv_file))


lema_proj   = crs(qudoba)
wgs         = "+init=EPSG:4326"
cols  = c("#FFFFFF","#92C5DE","#AFE1AF","#228B22","#808000")
gg_device  = c("png")     # Output devices to use (Check ggsave for acceptable formats)
gg_depth   = 600          # Plot resolution (dpi)
gg_ptsz    = 18           # Font size
gg_ptszl   = 26
gg_width   = 17.5         # Plot width (units below)
gg_widthn  = 14.5
gg_height  = 8.5          # Plot height (units below)
gg_units   = "in"         # Units for plot size
gg_screen  = TRUE         # Show plots on screen as well?
gg_tfmt    = "%y"         # Format for time strings in the time series %Y: 4 DIGITS YEAR %y 2 digits year
ndevice    = 1


quspba = raster::stack(qugaba,quchba,qudoba,quenba,quagba,qukeba,quloba,qumuba,quwiba)
sumba  = calc(quspba,sum)
sumba = sumba/100 #conver to m2/ha
#writeRaster(sumba,filename=file.path(main_lema,"CAoakba-sum.tif"),format="GTiff",overwrite=TRUE)
sumba = raster(file.path(main_lema,"CAoakba-sum-squaremeterPerHa.tif"))
ball  = ball/100 #conver to m2/ha
bacv  = bacv/100


prcnt1 = sumba/ball


ball_1 = ball
ball_1[ball_1<1] <- NA
prcnt2 = sumba/ball_1

sumba_1 = sumba
sumba_1[sumba_1<2.5] <- NA
prcnt3 = sumba_1/ball_1


oakba_prcnt = prcnt1


### resample to a raster similar to WRF domain
wrf_extnt = raster::extent(-130.2749,-108.9862,28.62024,45.71207)
wrf_poly  = as(wrf_extnt,"SpatialPolygons")
sp::proj4string(wrf_poly) = wgs
wrf_lema  = spTransform(wrf_poly,crs(lema_proj))
prcnt_wrf = crop(oakba_prcnt,extent(wrf_lema))
prcnt_wrf = mask(prcnt_wrf, wrf_lema)
prcnt_150 = aggregate(prcnt_wrf, fact=5, fun=mean)
prcnt_150[prcnt_150>1] <- NA
prcnt_wgs = projectRaster(prcnt_150, crs=wgs)



bacv_wrf  = crop(bacv, extent(wrf_lema))
bacv_wrf  = mask(bacv_wrf, wrf_lema)
bacv_150  = aggregate(bacv_wrf, fact=5,fun=mean)
bacv_wgs  = projectRaster(bacv_150, crs=wgs)

wrf_rs     = raster(wrf_extnt,ncols=147,nrows=151,crs=wgs)
prcnt_rspl = resample(prcnt_wgs, wrf_rs, method="bilinear")
bacv_rspl  = resample(bacv_wgs, wrf_rs, method="bilinear")

writeRaster(prcnt_rspl,filename=file.path(main_lema,"CAoakba-orcnt-WRF-All1Oak2p5.tif"),format="GTiff",overwrite=TRUE)


## resample to the model domain

bapnt_df = as.data.frame(prcnt_rspl,xy=TRUE)
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

wrf_qry        = wrf_df %>% dplyr::select(x_vec,y_vec,cellid,tobt_vec)

#search for the nearest neighbor in wrf for each corresponding lema cell
nn_pt        = RANN::nn2(bapnt_df[,1:2],wrf_qry[,1:2],k=1) 
wrf_qry$id   = as.vector(nn_pt$nn.idx) 
wrf_qry$dist = as.vector(nn_pt$nn.dists)

wrf_qry = wrf_qry                                     %>% 
          mutate(oakBA     = bapnt_df$layer[id])      %>% 
          rename(lon=x_vec,
                 lat=y_vec)                           %>% 
          dplyr::select(lon,lat,oakBA)

bacv_df      = as.data.frame(bacv_rspl,xy=TRUE)
cv_nn        = RANN::nn2(bacv_df[,1:2],wrf_qry[,1:2],k=1)
wrf_qry$id   = as.vector(cv_nn$nn.idx) 
wrf_qry$dist = as.vector(cv_nn$nn.dists)

wrf_qry = wrf_qry                                           %>% 
          mutate(cvBA     = bacv_df$ba_ge_3_2017_cv[id])    %>% 
          dplyr::select(-c(id,dist))


data.table::fwrite(wrf_qry,file.path(main_lema,"WRF_oakBasalAreaPrcnt_RatioLT1_CVmasked.csv"))

fil_var   = matrix(wrf_qry$oakBA, nrow=147,ncol=151)
x_arr     = matrix(wrf_qry$lon,nrow=147,ncol=151)
y_arr     = matrix(wrf_qry$lat,nrow=147,ncol=151)
wrf_star  = st_as_stars(fil_var)
wrf_star  = st_as_stars(wrf_star, curvilinear=list(X1=x_arr,X2=y_arr), crs=wgs)
wrf_sf    = st_as_sf(wrf_star,as_points=FALSE,na.rm=FALSE)
ca_co     = USAboundaries::us_counties(resolution = "high", states = "CA")

ba_plot = ggplot() + geom_sf(data=wrf_sf, aes(fill=A1),color=NA)+
  coord_sf(crs=wgs) + 
  scale_fill_gradientn(colours=cols, guide="colorbar", na.value="white") +
  geom_sf(data = ca_co, color = alpha("black", alpha=1),lwd=0.1,fill=NA) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) + 
  #labs(x=NULL,y=NULL,title=expression("LEMMA Oak BasalArea("~m^2~Ha~")")) +
  labs(x=NULL,title="LEMMA Oak BasalArea ratio")+
  theme_grey( base_size = gg_ptsz*0.8, base_family = "Helvetica",base_line_size= 0.5,base_rect_size =0.5) +
  theme(   legend.position  = "right"
           , legend.title     = element_blank()
           , legend.text      = element_text(size=gg_ptsz*0.8)
           , plot.title       = element_text(size=gg_ptsz*0.8, vjust=1)
           , panel.background = element_blank()
           , panel.border     = element_blank()
           , axis.text.x      = element_blank()
           , axis.text.y      = element_blank()
           , axis.ticks=element_blank()
           ,plot.margin = unit(c(0,0,0,0), "cm")
  )

for (d in sequence(ndevice)){
  h_output = paste0("lemma-oak-basalarea-0.7.",gg_device[d])
  dummy    = ggsave( filename = h_output
                     , plot      = ba_plot
                     , device    = gg_device[d]
                     , path      = main_lema
                     , width     = gg_widthn*0.6
                     , height    = gg_widthn*0.48
                     , units     = gg_units
                     , dpi       = gg_depth
  )
}


## canopy cover resample
cancov_wrf     = crop(cancov,extent(wrf_lema))
cancov_wrf     = mask(cancov_wrf, wrf_lema)
cancov150_wrf  = aggregate(cancov_wrf, fact=5, fun=mean)
cancov_wgs     = projectRaster(cancov150_wrf, crs=wgs)
cancov_mean    = resample(cancov_wgs, wrf_rs, method="bilinear")
## convert to raw units m2/ha
cancov_mean = cancov_mean/100
writeRaster(cancov_mean,filename=file.path(main_lema,"CA-HardwoodCanopycov-WRF.tif"),format="GTiff",
            overwrite=TRUE)

## resample to the model domain
cancov_df    = as.data.frame(cancov_mean,xy=TRUE)
wrf_qry      = wrf_df %>% dplyr::select(x_vec,y_vec,cellid)
nn_pt        = RANN::nn2(cancov_df[,1:2],wrf_qry[,1:2],k=1) 
wrf_qry$id   = as.vector(nn_pt$nn.idx) 
wrf_qry$dist = as.vector(nn_pt$nn.dists)

wrf_qry = wrf_qry                                               %>% 
          mutate(cancov    = cancov_df$cancov_hdw_2017[id])     %>% 
          rename(lon=x_vec,
                 lat=y_vec)                                     %>% 
  dplyr::select(lon,lat,cancov)
fwrite(wrf_qry,file.path(main_lema,"WRF_CA-HarwoodCanopyCover.csv"))

fil_var   = matrix(wrf_qry$cancov, nrow=147,ncol=151)
x_arr     = matrix(wrf_qry$lon,nrow=147,ncol=151)
y_arr     = matrix(wrf_qry$lat,nrow=147,ncol=151)
wrf_star  = st_as_stars(fil_var)
wrf_star  = st_as_stars(wrf_star, curvilinear=list(X1=x_arr,X2=y_arr), crs=wgs)
wrf_sf    = st_as_sf(wrf_star,as_points=FALSE,na.rm=FALSE)
ca_co     = USAboundaries::us_counties(resolution = "high", states = "CA")

cov_plot = ggplot() + geom_sf(data=wrf_sf, aes(fill=A1),color=NA)+
  coord_sf(crs=wgs) + 
  scale_fill_gradientn(colours=cols, guide="colorbar", na.value="white") +
  geom_sf(data = ca_co, color = alpha("black", alpha=1),lwd=0.1,fill=NA) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) + 
  labs(x=NULL,y=NULL,title="LEMMA Hardwood CanopyCover(%)") +
  theme_grey( base_size = gg_ptsz*0.8, base_family = "Helvetica",base_line_size= 0.5,base_rect_size =0.5) +
  theme(   legend.position  = "right"
           , legend.title     = element_blank()
           , legend.text      = element_text(size=gg_ptsz*0.8)
           , plot.title       = element_text(size=gg_ptsz*0.8, vjust=1)
           , panel.background = element_blank()
           , panel.border     = element_blank()
           , axis.text.x      = element_blank()
           , axis.text.y      = element_blank()
           , axis.ticks=element_blank()
           ,plot.margin = unit(c(0,0,0,0), "cm")
  )

for (d in sequence(ndevice)){
  h_output = paste0("lemma-hardwood-cancov.",gg_device[d])
  dummy    = ggsave( filename = h_output
                     , plot      = cov_plot
                     , device    = gg_device[d]
                     , path      = main_lema
                     , width     = gg_widthn*0.6
                     , height    = gg_widthn*0.48
                     , units     = gg_units
                     , dpi       = gg_depth
  )
}


## biomass resample
biomas_wrf     = crop(biomas,extent(wrf_lema))
biomas_wrf     = mask(biomas_wrf, wrf_lema)
biomas150_wrf  = aggregate(biomas_wrf, fact=5, fun=mean)
biomas_wgs     = projectRaster(biomas150_wrf, crs=wgs)
biomas_mean    = resample(biomas_wgs, wrf_rs, method="bilinear")
## convert from kg/ha to kgC/m2
biomas_mean    = biomas_mean*2*0.0001
writeRaster(biomas_mean,filename=file.path(main_lema,"CA-HardwoodBiomass-WRF.tif"),format="GTiff",
            overwrite=TRUE)

## resample to the model domain
biomas_df    = as.data.frame(biomas_mean,xy=TRUE)
wrf_qry      = wrf_df %>% dplyr::select(x_vec,y_vec,cellid)
nn_pt        = RANN::nn2(cancov_df[,1:2],wrf_qry[,1:2],k=1) 
wrf_qry$id   = as.vector(nn_pt$nn.idx) 
wrf_qry$dist = as.vector(nn_pt$nn.dists)

wrf_qry = wrf_qry                                       %>% 
  mutate(biomass    = biomas_df$bphh_ge_3_crm_2017[id]) %>% 
  rename(lon=x_vec,
         lat=y_vec)                                     %>% 
  dplyr::select(lon,lat,biomass)
fwrite(wrf_qry,file.path(main_lema,"WRF_CA-HarwoodBiomass.csv"))

fil_var   = matrix(wrf_qry$biomass, nrow=147,ncol=151)
x_arr     = matrix(wrf_qry$lon,nrow=147,ncol=151)
y_arr     = matrix(wrf_qry$lat,nrow=147,ncol=151)
wrf_star  = st_as_stars(fil_var)
wrf_star  = st_as_stars(wrf_star, curvilinear=list(X1=x_arr,X2=y_arr), crs=wgs)
wrf_sf    = st_as_sf(wrf_star,as_points=FALSE,na.rm=FALSE)
ca_co     = USAboundaries::us_counties(resolution = "high", states = "CA")

ms_plot = ggplot() + geom_sf(data=wrf_sf, aes(fill=A1),color=NA)+
  coord_sf(crs=wgs) + 
  scale_fill_gradientn(colours=cols, guide="colorbar", na.value="white") +
  geom_sf(data = ca_co, color = alpha("black", alpha=1),lwd=0.1,fill=NA) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) + 
  labs(x=NULL,y=NULL,title=expression("LEMMA Hardwood Biomass("~kgC~m^-2~")")) +
  theme_grey( base_size = gg_ptsz*0.8, base_family = "Helvetica",base_line_size= 0.5,base_rect_size =0.5) +
  theme(   legend.position  = "right"
          , legend.title     = element_blank()
          , legend.text      = element_text(size=gg_ptsz*0.8)
          , plot.title       = element_text(size=gg_ptsz*0.8, vjust=1)
          , panel.background = element_blank()
          , panel.border     = element_blank()
          , axis.text.x      = element_blank()
          , axis.text.y      = element_blank()
          , axis.ticks=element_blank()
          ,plot.margin = unit(c(0,0,0,0), "cm")
  )

for (d in sequence(ndevice)){
  h_output = paste0("lemma-hardwood-biomass.",gg_device[d])
  dummy    = ggsave( filename = h_output
                  , plot      = ms_plot
                  , device    = gg_device[d]
                  , path      = main_lema
                  , width     = gg_widthn*0.6
                  , height    = gg_widthn*0.48
                  , units     = gg_units
                  , dpi       = gg_depth
  )
}

#### develope a strip domain along the precipitation gradient 
### which has canopy cover varying 10% - 40%
ba = data.table::fread("~/Google Drive/My Drive/CA-OakWoodland/benchmark/WRF_oakBasalArea.csv")
strip = ba                                                                              %>% 
        mutate(mask1 = ifelse(lat>=37.5 & lat<=37.8 & lon>=-120.95 & lon<=-119.73,1,0)) %>% 
        mutate(mask  = ifelse(mask1==1 & oakBA > 0, 1, 0) )                             %>% 
        select(lon,lat,mask,oakBA)                                                      %>% 
        mutate(oakBA=ifelse(mask==1,oakBA,0))

fil_var   = matrix(strip$oakBA, nrow=147,ncol=151)
x_arr     = matrix(strip$lon,nrow=147,ncol=151)
y_arr     = matrix(strip$lat,nrow=147,ncol=151)
wrf_star  = st_as_stars(fil_var)
wrf_star  = st_as_stars(wrf_star, curvilinear=list(X1=x_arr,X2=y_arr), crs=wgs)
wrf_sf    = st_as_sf(wrf_star,as_points=FALSE,na.rm=FALSE)
ca_co     = USAboundaries::us_counties(resolution = "high", states = "CA")

stp_plot = ggplot() + geom_sf(data=wrf_sf, aes(fill=A1),color=NA)+
  coord_sf(crs=wgs) + 
  scale_fill_gradientn(colours=cols, guide="colorbar", na.value="white") +
  geom_sf(data = ca_co, color = alpha("black", alpha=1),lwd=0.1,fill=NA) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) + 
  labs(x=NULL,y=NULL,title=expression("LEMMA Oak BasalArea ("~m^2~Ha~")")) +
  theme_grey( base_size = gg_ptsz*0.8, base_family = "Helvetica",base_line_size= 0.5,base_rect_size =0.5) +
  theme(   legend.position  = "right"
           , legend.title     = element_blank()
           , legend.text      = element_text(size=gg_ptsz*0.8)
           , plot.title       = element_text(size=gg_ptsz*0.8, vjust=1)
           , panel.background = element_blank()
           , panel.border     = element_blank()
           , axis.text.x      = element_blank()
           , axis.text.y      = element_blank()
           , axis.ticks=element_blank()
           ,plot.margin = unit(c(0,0,0,0), "cm")
  )


for (d in sequence(ndevice)){
  h_output = paste0("oakwoodland_domain.",gg_device[d])
  dummy    = ggsave( filename = h_output
                     , plot      = stp_plot
                     , device    = gg_device[d]
                     , path      = main_lema
                     , width     = gg_widthn*0.6
                     , height    = gg_widthn*0.48
                     , units     = gg_units
                     , dpi       = gg_depth
  )
}


## create domain mask
author_name  = "Xiulin Gao"
author_email = "xiulingao@lbl.gov"
undef        = -9999.0000

land_mask  <- array(strip$mask1,dim=c(147,151))
land_mkdif <- land_mask


## create new nc file as land mask
xx  = ncdim_def( name="lon"   ,units="",vals= sequence(147)  ,create_dimvar=FALSE)
yy  = ncdim_def( name="lat"   ,units="",vals= sequence(151)  ,create_dimvar=FALSE)
nc_xy  = list   (xx,yy)
xy     = c(147,151)
file_name = file.path("Google Drive/My Drive/CA-OakWoodland/strip-mask.nc")
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

att_template = list( title            = "To be replaced when looping through months"
                     , date_created   = paste0(as.character(now(tzone="UTC")), "UTC")
                     , source_code    = "mask-wrf-domain.R"
                     , code_notes     = "land mask for WRF domain created a strip run in CA oak woodland"
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


nc_title   = "Land mask for a strip run across precipitation gradient in CA oak woodland"
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



### oak woodland domain using the oak basal area%
bapnt = data.table::fread(file.path(main_lema,"WRF_oakBasalAreaPrcnt_RatioLT1_CVmasked.csv"))
domain = bapnt                                            %>%
         mutate(mask1 = ifelse(oakBA<0.5, 0, 1),
                mask2 = ifelse(cvBA<0.6 & !is.na(cvBA),0,1))    %>% 
         mutate(mask = ifelse(mask1 & mask2, 1, 0))
domain = domain                                                %>% 
         mutate(mask= ifelse(lat > 41, 0,mask),
                mask = ifelse(lat < 34.3 & lat > 33 & lon > -121 & lon < -119.6 , 0, mask),
                mask = ifelse(lon > -116, 0, mask),
                mask = ifelse(lat > 39.5 & lat < 41 & lon < -119 & lon > -121, 0, mask),
                mask = ifelse(lat > 35 & lat < 37 & lon > -117.5 & lon < -115, 0, mask))
domain = domain %>% mutate(mask = ifelse(is.na(mask), 0, mask))


fil_var   = matrix(as.factor(domain$mask), nrow=147,ncol=151)
x_arr     = matrix(domain$lon,nrow=147,ncol=151)
y_arr     = matrix(domain$lat,nrow=147,ncol=151)
wrf_star  = st_as_stars(fil_var)
wrf_star  = st_as_stars(wrf_star, curvilinear=list(X1=x_arr,X2=y_arr), crs=wgs)
wrf_sf    = st_as_sf(wrf_star,as_points=FALSE,na.rm=FALSE)
ca_co     = USAboundaries::us_counties(resolution = "high", states = "CA")

stp_plot = ggplot() + geom_sf(data=wrf_sf, aes(fill=A1),color=NA,lwd=0)+
  coord_sf(crs=wgs) + #geom_point(aes(x=-122.69,y=38.57), colour=alpha("black",1), size=1)+
  scale_fill_manual(values=c("1" = "lightgreen","0"="white"),na.value="white") +
  geom_sf(data = ca_co, color = alpha("black", alpha=1),lwd=0.1,fill=NA) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) + 
  labs(x=NULL,y=NULL,title=NULL) +
  theme_grey( base_size = gg_ptsz*0.8, base_family = "Helvetica",base_line_size= 0.5,base_rect_size =0.5) +
  theme(     legend.position  = "none"
           #, legend.title     = element_blank()
           #, legend.text      = element_text(size=gg_ptsz*0.8)
           #, plot.title       = element_text(size=gg_ptsz*0.8, vjust=1)
           , panel.background = element_blank()
           , panel.border     = element_blank()
           , rect = element_blank()
           , axis.text.x      = element_blank()
           , axis.text.y      = element_blank()
           , axis.ticks=element_blank()
           ,plot.margin = unit(c(0,0,0,0), "cm")
  )

stp_plot


for (d in sequence(ndevice)){
  h_output = paste0("oakwoodland_domain.",gg_device[d])
  dummy    = ggsave( filename = h_output
                     , plot      = stp_plot
                     , device    = gg_device[d]
                     , path      = main_lema
                     , width     = gg_widthn*0.6
                     , height    = gg_widthn*0.48
                     , units     = gg_units
                     , dpi       = gg_depth
  )
}


land_mask  <- array(domain$mask,dim=c(147,151))
land_mkdif <- land_mask


## create new nc file as land mask
xx  = ncdim_def( name="lon"   ,units="",vals= sequence(147)  ,create_dimvar=FALSE)
yy  = ncdim_def( name="lat"   ,units="",vals= sequence(151)  ,create_dimvar=FALSE)
nc_xy  = list   (xx,yy)
xy     = c(147,151)
file_name = file.path("Google Drive/My Drive/CA-OakWoodland/oakwoodland-domain-mask-BAratio-oakBALT1GTP5.nc")
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

att_template = list( title            = "To be replaced when looping through months"
                     , date_created   = paste0(as.character(now(tzone="UTC")), "UTC")
                     , source_code    = "mask-wrf-domain.R"
                     , code_notes     = "land mask for CA oak woodland created using LEMMA oak basal area to total basal area > 0.5)"
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


nc_title   = "Land mask for CA oak woodland using LEMMA oak basal to total basal area >0.5"
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
