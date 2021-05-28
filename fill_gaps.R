# Stephanie.Clay@dfo-mpo.gc.ca

# rm(list=ls())

library(metR)
library(dplyr)
library(data.table)
library(fst)
library(oceancolouR)
library(lubridate)
library(geodist) # to create a large continuous artificial gap in the image
source("ImputeEOF2.R")
source("load_and_subset.R")


region <- "PANCAN"

# For MODIS/SeaWiFS/VIIRS-SNPP: CHL_OCX, CHL_POLY4, CHL_GSM_GS, PAR
variable <- "CHL_OCX"

# MODIS, VIIRS-SNPP, SeaWiFS
sensor <- "MODIS"

# year to fill
year <- 2015

path <- "/mnt/data3/claysa"
output_path <- "method_testing_output"


#**************
# QUALITY CONTROL OPTIONS

# remove days with less than this percent coverage
low_percov <- 5

# set boundaries of valid chla to use
# (note that if you're filling gaps in logged chla, values <= 0 will be set to NA automatically)
low_chla <- 0
high_chla <- 64

#**************

# for 8day, how many years to use on either side of the target year?
# if filling gaps for weekly data, you can load a few years at a time
# if filling gaps in daily data, you can only load ~ 1 year
num_years <- 0 # default = 4 (for 8day filling)

# log the data before filling? (for CHL)
fill_log <- TRUE

# fill gaps in 8day or daily images
composite <- "daily"



# # note that you can comment out fill_log and composite above, and run this code in the console instead:
# for (fill_log in c(TRUE, FALSE)) {
#     for (composite in c("8day", "daily")) {
#         source("fill_gaps.R")
#     }
# }




#*******************************************************************************
# MAIN CODE ####

output_years <- ifelse(num_years==0, year, paste0(range(plus_minus(year,num_years)), collapse="-"))
output_name <- paste0(output_path, "/filled_", region, "_", sensor, "_", variable, "_", output_years, "_", composite,
                      ifelse(fill_log, "_logged", ""), "_randomCVpts")

data("pancan_lats_4km")
data("pancan_lons_4km")
data("pancan_bins_4km")
data("nwa_lats_4km")
data("nwa_lons_4km")
data("nwa_bins_4km")
data("nep_lats_4km")
data("nep_lons_4km")
data("nep_bins_4km")

# get 4km-resolution lats/lons in vector form
if (region=="PANCAN") {
    input_lats <- pancan_lats_4km
    input_lons <- pancan_lons_4km
    input_bins <- pancan_bins_4km
    output_bins <- nwa_bins_4km
    output_lats <- nwa_lats_4km
    output_lons <- nwa_lons_4km
} else if (region=="NWA") {
    input_lats <- nwa_lats_4km
    input_lons <- nwa_lons_4km
    input_bins <- nwa_bins_4km
    output_bins <- nwa_bins_4km
    output_lats <- nwa_lats_4km
    output_lons <- nwa_lons_4km
} else if (region=="NEP") {
    input_lats <- nep_lats_4km
    input_lons <- nep_lons_4km
    input_bins <- nep_bins_4km
    output_bins <- nep_bins_4km
    output_lats <- nep_lats_4km
    output_lons <- nep_lons_4km
}

# collect attributes (number of eof, rmse...)
attr_filename <- paste0(output_name, "_attr.rda")
if (file.exists(file.path(path, attr_filename))) {
    load(file.path(path, attr_filename))
} else {
    attr_df <- data.frame(sensor=character(),
                          variable=character(),
                          region=character(),
                          year=integer(),
                          composite=character(),
                          fill_log=logical(),
                          eof=numeric(),
                          rmse=numeric(),
                          num_LE_zero=integer(),
                          processing_time=character(),
                          time_to_process=numeric(),
                          stringsAsFactors=FALSE)
}


ptm <- Sys.time()

# get input files
file_list <- sort(list.files(file.path(path, sensor, variable, region, "annual_fst")))
filelist_years <- as.numeric(substr(sapply(sapply(file_list, strsplit, split="_"), "[[", 5),1,4))
file_list <- file_list[filelist_years %in% plus_minus(year,num_years)]

print(file_list)

# load and format input files
df <- lapply(file_list,
             load_and_subset,
             path=file.path(path, sensor, variable, region, "annual_fst"),
             input_bins=input_bins,
             output_lats=output_lats,
             output_lons=output_lons,
             output_bins=output_bins,
             low_percov=low_percov,
             temporal_resolution=composite)
df <- do.call(dplyr::bind_rows, df)

# set values outside acceptable boundaries to NA
df[is.finite(df$var) & (df$var < low_chla | df$var > high_chla), "var"] <- NA

if (fill_log) {
    good_ind <- is.finite(df$var) & df$var > 0
    df[!good_ind,"var"] <- NA
    df[good_ind,"var"] <- log10(df[good_ind,"var"])
}

# get existing CV pixels from the main year
if (num_years > 0) {
    CV_pixels <- read_fst(paste0(output_path, "/filled_", region, "_", sensor, "_", variable, "_", year, "_", composite,
                                 ifelse(fill_log, "_logged", ""), "_randomCVpts.fst")) %>%
        dplyr::filter(floor(time)==year)
    CV_pixels <- CV_pixels$validation_pixels
    all_years <- sort(unique(floor(df$time)))
    available_time <- lapply(all_years, function(y) as.integer(round(sort(unique(((df %>% dplyr::filter(floor(time)==y))$time-y)*46 + 1)))))
    filled_data <- ImputeEOF2(formula = "var ~ lon + lat | time",
                              data = df,
                              num_rows = length(output_bins),
                              tol = ifelse(fill_log, 0.001, 0.01),
                              CV_pixels = CV_pixels,
                              all_years=all_years,
                              available_time=available_time)
} else {
    filled_data <- ImputeEOF2(formula = "var ~ lat + lon | time",
                              data = df,
                              num_rows = length(output_bins),
                              tol = ifelse(fill_log, 0.001, 0.01),
                              all_years=year)
}

df$var_imputed <- filled_data$X_filled
df$validation_pixels <- filled_data$validation_pixels

if (fill_log) {
    df <- df %>% dplyr::mutate(var = 10^var,
                               var_imputed = 10^var_imputed,
                               validation_pixels = 10^validation_pixels)
}

# get the number of filled values that are <= 0, for reporting
num_LE_zero <- sum(df$var_imputed <= 0, na.rm=TRUE)

attr_df <- dplyr::bind_rows(attr_df,
                            data.frame(sensor=as.character(sensor),
                                       variable=variable,
                                       region=region,
                                       year=year,
                                       composite=composite,
                                       fill_log=fill_log,
                                       eof=attr(df$var_imputed, "eof"),
                                       rmse=attr(df$var_imputed, "rmse"),
                                       num_LE_zero=num_LE_zero,
                                       processing_time=format(Sys.time(), "%Y%m%d_%H%M%S"),
                                       time_to_process=as.numeric(Sys.time()-ptm, units="secs"),
                                       stringsAsFactors = FALSE))

write_fst(df %>% dplyr::select(bin, time, var, var_imputed, validation_pixels),
          path=paste0(output_name, ".fst"),
          compress=100)

# clear memory
rm("df", "filled_data")
gc()

# save attributes (eof number and rmse) from each call of ImputeEOF
save(attr_df, file=attr_filename, compress=TRUE)

