#####################################################################################################################
# This script produces:
# A) WSC MONTHLY AND ANNUAL DATA (1981–2023) spatially explicit FOR ROI: lon(-48, 8), lat(40, 75) ----
            # 1) df_monthly_roi      – WSC monthly avg. (1981–2023) over ROI
            # 2) df_annual_roi       – WSC annual avg. (1981–2023) over ROI
# B) WSClag1 INDICES (1981–2023) used in regression models (averaged over BOX: lon(-17, -12), lat(52, 55)) ----
            # 1) WSClag1_JanDec_DF   – WSC lag-1 annual (Jan–Dec) avg. over BOX for 1981-2023
            # 2) WSClag1_MAM_DF      – WSC lag-1 spawning season (Mar–May) avg. over BOX for 1981-2023
#####################################################################################################################




# Required packages ----
library(ncdf4)       
library(dplyr)     
library(tidyr)       
library(lubridate)   
library(openxlsx)    
library(reshape2)    


# Set working directory ----
#setwd("")


##################################################################################################----
########## PREPARE WIND STRESS DATA: read NetCDF files and compute monthly wind stress curl (1948–2024) ----
##################################################################################################----
# NetCDF files (uwnd, vwnd): wind stress monthly avg, 10m resolution 
# uwnd = Eastward wind component; vwnd = Northward wind component
# Source: NCAR/NCEP Reanalysis project
# Note: data downloaded in October 2024, so contains data available for January 1948 - September 2024
uwnd <- nc_open("Data/windstress_NCARNCEPreanalysis_20241002/uwnd10m.mon.mean.nc")
vwnd <- nc_open("Data/windstress_NCARNCEPreanalysis_20241002/vwnd10m.mon.mean.nc")

# Get variables from NetCDF files
lat <- ncvar_get(uwnd, 'lat')
lon <- ncvar_get(uwnd, 'lon')
time <- ncvar_get(uwnd, 'time')
uflx <- ncvar_get(uwnd, 'uwnd')
vflx <- ncvar_get(vwnd, 'vwnd')

# Close the NetCDF files
nc_close(uwnd)
nc_close(vwnd)

# Permute dimensions (order dimensions as [time, lat, lon])
uflx <- aperm(uflx, c(3, 2, 1))
vflx <- aperm(vflx, c(3, 2, 1))
# Adjust longitudes and wind arrays to the required format
lon <- c(lon[(length(lon)/2 + 1):length(lon)] - 360, lon[1:(length(lon)/2)])
uflx <- uflx[, , c((length(lon)/2 + 1):length(lon), 1:(length(lon)/2))]
vflx <- vflx[, , c((length(lon)/2 + 1):length(lon), 1:(length(lon)/2))]

# Grid setup
Nt <- dim(uflx)[1]; Ny <- dim(uflx)[2]; Nx <- dim(uflx)[3]
LON <- matrix(rep(lon, each = Ny), nrow = Ny)
LAT <- matrix(rep(lat, times = Nx), ncol = Nx)
dx <- cos(LAT * pi / 180) * 111111
dx <- dx[, 1:(Nx - 1)]
dy <- matrix(111111, nrow = Ny - 1, ncol = Nx)

# Compute curl of wind stress: ∂v/∂x − ∂u/∂y
# Where u = eastward wind stress, v = northward wind stress
#       dvdx ≈ (v[x+1] − v[x]) / dx ; dudy ≈ (u[y+1] − u[y]) / dy
# Curl = dvdx - dudy
cav <- array(0, dim = c(Nt, Ny, Nx))
for (t in 1:Nt) {
  u <- -uflx[t, , ]; v <- -vflx[t, , ]
  dvdx <- (v[, 2:Nx] - v[, 1:(Nx - 1)]) / dx
  dudy <- -(u[2:Ny, ] - u[1:(Ny - 1), ]) / dy
  dvdx <- cbind(dvdx, dvdx[, 1])
  dudy <- rbind(dudy, dudy[1, ])
  cav[t, , ] <- dvdx - dudy
}
# WSC monthly dataframe
wsc_monthly_CC_ <- list(lat = matrix(lat, ncol = 1),
                        lon = matrix(lon, ncol = 1),
                        time = matrix(time, ncol = 1),
                        wsc = cav)

# Convert time to date format
time_units <- "hours since 1800-01-01 00:00:0.0"
days_since_ref <- time / 24
reference_date <- as.Date("1800-01-01")
tt <- reference_date + days_since_ref


##################################################################################################----
########## A) EXTRACT WSC MONTHLY AND ANNUAL DATA (1981–2023) FOR ROI: lon(-48, 8), lat(40, 75) ----
##################################################################################################----

# Rearrange arrays and adjust sign
wsc_monthly_CC <- aperm(-wsc_monthly_CC_$wsc, c(3, 2, 1))[, 94:1, ]

lat_windstress <- rev(wsc_monthly_CC_$lat)
lon_windstress <- wsc_monthly_CC_$lon

date_windstress <- as.POSIXct(wsc_monthly_CC_$time * 60 * 60, origin = "1800-01-01")

# Define ROI
lons_ind_map <- which(lon_windstress >= -48 & lon_windstress <= 8)
lats_ind_map <- which(lat_windstress > 40 & lat_windstress < 75)

# Extract WSC data for ROI
roi_data <- wsc_monthly_CC[lons_ind_map, lats_ind_map, ]
# Assign dimnames to the array 
dimnames(roi_data) <- list(
  longitude = as.character(lon_windstress[lons_ind_map]),
  latitude = as.character(lat_windstress[lats_ind_map]),
  date = as.character(date_windstress[1:dim(roi_data)[3]]))
# Convert to data frame
roi_df <- as.data.frame.table(roi_data)
names(roi_df)[4] <- "wsc"
# Convert date string to POSIXct
roi_df$date <- as.POSIXct(roi_df$date, tz = "UTC")
# Add year and month columns
roi_df$year <- year(roi_df$date)
roi_df$month <- month(roi_df$date)
# Add lon, lat 
roi_df$longitude <- as.numeric(as.character(roi_df$longitude))
roi_df$latitude <- as.numeric(as.character(roi_df$latitude))
# Keep only 1981-2023 data
roi_df <- roi_df %>% filter(year >= 1981 & year <= 2023)

# --- 1) WSC monthly avg. for 1981-2023 over ROI ----
df_monthly_roi <- roi_df

# --- 2) WSC annual avg. for 1981-2023 over ROI ----
df_annual_roi <- roi_df %>%
  group_by(year, longitude, latitude) %>%
  summarise(wsc = mean(wsc, na.rm = TRUE)) %>%
  ungroup()



##################################################################################################----
########## B) EXTRACT WSClag1 INDICES (1981–2023) used in regression models (BOX: lon(-17, -12), lat(52, 55))  ----
##################################################################################################----

# Define BOX area
lons_ind_box <- which(lon_windstress >= -17 & lon_windstress <= -12)
lats_ind_box <- which(lat_windstress > 52 & lat_windstress < 55)
lons_box <- lon_windstress[lons_ind_box]
lats_box <- lat_windstress[lats_ind_box]

monthly_box_data <- wsc_monthly_CC[lons_ind_box, lats_ind_box, ]
dimnames(monthly_box_data) <- list(longitude = as.character(lons_box),
                                   latitude = as.character(lats_box),
                                   date = as.character(date_windstress[1:dim(monthly_box_data)[3]]))

# Convert to data frame
box_df <- as.data.frame.table(monthly_box_data)
names(box_df)[4] <- "wsc"
box_df$date <- as.POSIXct(box_df$date, tz = "UTC")
box_df$year <- year(box_df$date)
box_df$month <- month(box_df$date)

# Monthly mean WSC over the BOX
WSC_monthly_BOX_DF <- box_df %>%
  group_by(year, month) %>%
  summarise(monthly_mean_wsc = mean(wsc, na.rm = TRUE)) %>%
  ungroup()


# --- 1) WSClag1 annual (Jan–Dec) average over BOX (lag-1: WSC Jan–Dec in 1980 → WSClag1 in 1981) ----
WSClag1_JanDec_DF <- WSC_monthly_BOX_DF %>%
  filter(year >= 1980, year <= 2022, month %in% 1:12) %>%
  group_by(year) %>%
  summarise(mean_wsc = mean(monthly_mean_wsc, na.rm = TRUE)) %>%
  mutate(year = year + 1) %>%
  filter(year >= 1981, year <= 2023) %>%
  rename(WSClag1_JanDec = mean_wsc)

# --- 2) WSClag1 spawning season (Mar–May) average over BOX (lag-1: WSC Mar–May in 1980 → WSClag1 in 1981) ----
WSClag1_MAM_DF <- WSC_monthly_BOX_DF %>%
  filter(year >= 1980, year <= 2022, month %in% 3:5) %>%
  group_by(year) %>%
  summarise(WSClag1_MAM = mean(monthly_mean_wsc, na.rm = TRUE)) %>%
  mutate(year = year + 1) %>%
  filter(year >= 1981, year <= 2023)


