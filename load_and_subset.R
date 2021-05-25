# function to load an annual fst file, remove days/weeks with low percent coverage,
# and subset it to the desired region if necessary
load_and_subset <- function(f, path, input_bins,
                            output_lats, output_lons, output_bins=input_bins,
                            low_percov=0, temporal_resolution="daily") {
    
    fsplit <- as.numeric(substr(strsplit(f, "_")[[1]], 1, 4))
    year <- as.numeric(fsplit[is.finite(fsplit)])
    
    sat_data <- read_fst(file.path(path, f))
    
    # read input data and reduce to output bins if necessary
    if (!identical(input_bins,output_bins)) {
        num_input_days <- nrow(sat_data)/length(input_bins)
        sat_data <- sat_data %>%
            dplyr::mutate(bin = rep(input_bins, num_input_days)) %>%
            dplyr::filter(bin %in% output_bins) %>%
            dplyr::select(var)
    }

    sat_mat <- matrix(sat_data$var, nrow=length(output_bins))
    
    if (temporal_resolution=="daily") {
        full_num_composites <- ifelse(leap_year(year),366,365)
    } else if (temporal_resolution=="8day") {
        full_num_composites <- 46
        sat_mat <- convert_daily_grid(sat_mat, composite="8day")
    }
    
    num_composites <- ncol(sat_mat)
    
    good_percov_ind <- 100 * colSums(is.finite(sat_mat)) / length(output_bins) >= low_percov
    
    df <- data.frame(bin = rep(output_bins, sum(good_percov_ind)),
                     lat = rep(output_lats, sum(good_percov_ind)),
                     lon = rep(output_lons, sum(good_percov_ind)),
                     time = year + (rep((1:num_composites)[good_percov_ind], each=length(output_bins)) - 1)/full_num_composites,
                     var = as.numeric(sat_mat[,good_percov_ind]),
                     stringsAsFactors = FALSE)
    
    return(df)
    
}
