# function to load an annual fst file, remove days/weeks with low percent coverage,
# and subset it to the desired region if necessary
load_and_subset <- function(f, path, input_region, input_bins, input_lats, input_lons,
                            output_bins=NULL, low_percov=0, spatial_resolution="4km",
                            temporal_resolution="daily") {
    
    fsplit <- as.numeric(substr(strsplit(f, "_")[[1]], 1, 4))
    year <- as.numeric(fsplit[is.finite(fsplit)])
    
    sat_data <- read_fst(file.path(path, f))
    sat_mat <- matrix(sat_data$var, nrow=num_pix[[input_region]][[spatial_resolution]])
    
    if (temporal_resolution=="daily") {
        full_num_composites <- ifelse(leap_year(year),366,365)
    } else if (temporal_resolution=="8day") {
        full_num_composites <- 46
        sat_mat <- convert_daily_grid(sat_mat, composite="8day")
    }
    
    num_composites <- ncol(sat_mat)
    
    good_percov_ind <- 100 * colSums(is.finite(sat_mat)) / num_pix[[input_region]][[spatial_resolution]] >= low_percov
    df <- dplyr::left_join(data.frame(bin = output_bins, stringsAsFactors = FALSE),
                           data.frame(var = as.numeric(sat_mat[,good_percov_ind])) %>%
                               dplyr::mutate(bin = rep(input_bins, sum(good_percov_ind)),
                                             lat = rep(input_lats, sum(good_percov_ind)),
                                             lon = rep(input_lons, sum(good_percov_ind)),
                                             time = year + (rep((1:num_composites)[good_percov_ind], each=num_pix[[input_region]][[spatial_resolution]]) - 1)/full_num_composites),
                           by="bin") %>%
        as.data.table()
    
    return(df)
    
}
