
# make map comparing original and filled images
make_map <- function(df,title,log_raster=FALSE) {
    r1 <- var_to_rast(df %>% dplyr::select(bin, var))
    r2 <- var_to_rast(df %>% dplyr::select(bin, var_imputed))
    rextent <- raster::extent(c(xmn=lons[1], xmx=lons[2], ymn=lats[1], ymx=lats[2]))
    r1 <- crop(r1, rextent)
    r2 <- crop(r2, rextent)
    if (log_raster) {
        r1 <- log10(r1)
        r2 <- log10(r2)
    }
    r <- raster::stack(r1,r2)
    names(r) <- c("Original", "Filled")
    return(spplot(r, main=title, zlim=c(-1.5,2)) + latticeExtra::layer(sp::sp.polygons(wrld_simpl)))
}


# make plot of cross-validation points and their regression line
cv_plot <- function(cv_df,intercept,slope,tab,title) {
    breaks <- c(0.01,0.1,1,10,100)
    ggplot(cv_df, aes(x=var, y=validation_pixels)) +
        stat_binhex(bins=80) +
        scale_fill_gradientn(colours=c("#DDDDDD", "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), na.value="#000000FF") +
        scale_x_log10(limits=range(breaks),breaks=breaks,labels=breaks) +
        scale_y_log10(limits=range(breaks),breaks=breaks,labels=breaks) +
        geom_abline(intercept=0, slope=1) +
        geom_abline(intercept=intercept, slope=slope, colour="red", alpha=0.6, linetype="dashed") +
        labs(x="Real value", y="Filled value") +
        ggtitle(title) +
        theme_bw() +
        theme(legend.position="none") +
        annotation_custom(tab, xmin=-Inf, xmax=0.1, ymin=0.8, ymax=Inf)
}


# make table of regression stats for CV (cross-validation) regression and in situ/satellite regression
make_tab <- function(y,x) {
    tmp_lm <- lmodel2(y ~ x)
    stats_df <- data.frame(matrix(c(as.numeric(broom::glance(tmp_lm)[,c("r.squared", "p.value", "nobs")]),
                                    as.numeric(unlist(broom::tidy(tmp_lm)[5:6,"estimate"]))), ncol=1),
                           row.names=c("R^2", "p-value", "nobs", "intercept", "slope"),
                           stringsAsFactors = FALSE)
    tab <- tableGrob(round(stats_df,3), cols=NULL)
    return(list(stats_df=stats_df,tab=tab))
}


# make plot of in situ/satellite matchups
matchup_plot <- function(matchups) {
    lm_insitu <- make_tab(y=log10(matchups$sat_chla), x=log10(matchups$HPLC_CHLA))
    lm_real <- broom::tidy(lmodel2(log10(sat_chla) ~ log10(HPLC_CHLA), data=matchups %>% dplyr::filter(Satellite=="Real")))
    lm_filled <- broom::tidy(lmodel2(log10(sat_chla) ~ log10(HPLC_CHLA), data=matchups %>% dplyr::filter(Satellite=="Filled")))
    
    ggplot(matchups, aes(x=HPLC_CHLA, y=sat_chla, color=Satellite)) +
        geom_point() +
        theme_bw() +
        geom_abline(intercept=0, slope=1) +
        geom_abline(intercept=lm_insitu$stats_df[4,], slope=lm_insitu$stats_df[5,], linetype="dashed") +
        geom_abline(intercept=lm_filled$estimate[5], slope=lm_filled$estimate[6], color="#F8766D", linetype="dashed") +
        geom_abline(intercept=lm_real$estimate[5], slope=lm_real$estimate[6], color="#00B0F6", linetype="dashed") +
        scale_x_log10(limits=c(0.05,20)) +
        scale_y_log10(limits=c(0.05,20)) +
        labs(x="In situ chla", y="Satellite chla", color="") +
        annotation_custom(lm_insitu$tab, xmin=-Inf, xmax=0.06, ymin=0.5, ymax=Inf) +
        theme(legend.position=c(0.85,0.15))
}


# given a number of years used in the reconstruction, return some maps, plots, and stats to gauge the performance
time_series_comparison <- function(variable, region, fill_log, composite, num_years) {
    input_years <- ifelse(num_years==0, year, paste0(range(plus_minus(year,num_years)), collapse="-"))
    input_name <- paste0("output/filled_", region, "_", sensor, "_", variable, "_", input_years, "_")
    
    # load data
    file <- paste0(input_name, composite, ifelse(fill_log, "_logged", ""), "_randomCVpts")
    df <- as.data.table(read_fst(paste0(file, ".fst")))
    load(paste0(file, "_attr.rda"))
    
    # format
    num_weeks <- dim(df)[1]/num_pix$NWA$`4km`
    available_weeks <- as.integer(round(sort(unique((df$time-year)*46 + 1))))
    mat_var <- matrix(df$var, ncol=num_weeks)
    mat_var_imputed <- matrix(df$var_imputed, ncol=num_weeks)
    
    # make maps for each week
    map1 <- make_map(data.frame(bin=nwa_bins_4km,
                                var=mat_var[,available_weeks==weeks[1]],
                                var_imputed=mat_var_imputed[,available_weeks==weeks[1]],
                                stringsAsFactors = FALSE),
                     title=paste0("8day filled - week ",weeks[1]),
                     log_raster=TRUE)
    map2 <- make_map(data.frame(bin=nwa_bins_4km,
                                var=mat_var[,available_weeks==weeks[2]],
                                var_imputed=mat_var_imputed[,available_weeks==weeks[2]],
                                stringsAsFactors = FALSE),
                     title=paste0("8day filled - week ",weeks[2]),
                     log_raster=TRUE)
    
    # perform linear regression on all cross-validation pixels, and extract stats
    cv_df <- df %>% dplyr::filter(is.finite(validation_pixels))
    tab <- make_tab(y=log10(cv_df$validation_pixels), x=log10(cv_df$var))
    p <- cv_plot(cv_df,tab$stats_df[4,],tab$stats_df[5,],tab$tab,"CV regression - all pixels")
    
    # do the same with individual selected weeks to test
    week1_cv_df <- cv_df %>% dplyr::filter(time==(year + (weeks[1]-1)/46))
    week2_cv_df <- cv_df %>% dplyr::filter(time==(year + (weeks[2]-1)/46))
    if (fill_log) {
        week_rmses <- c(rmse(log10(week1_cv_df$var), log10(week1_cv_df$validation_pixels)),
                        rmse(log10(week2_cv_df$var), log10(week2_cv_df$validation_pixels)))
    } else {
        week_rmses <- c(rmse(week1_cv_df$var, week1_cv_df$validation_pixels),
                        rmse(week2_cv_df$var, week2_cv_df$validation_pixels))
    }
    tab <- make_tab(y=log10(week1_cv_df$validation_pixels), x=log10(week1_cv_df$var))
    p1 <- cv_plot(week1_cv_df,tab$stats_df[4,],tab$stats_df[5,],tab$tab,paste0("CV regression - week ", weeks[1]))
    tab <- make_tab(y=log10(week2_cv_df$validation_pixels), x=log10(week2_cv_df$var))
    p2 <- cv_plot(week2_cv_df,tab$stats_df[4,],tab$stats_df[5,],tab$tab,paste0("CV regression - week ", weeks[2]))
    
    # remove repetitive axis titles/scales
    p2 <- p2 + theme(axis.title.y=element_blank(), axis.text.y=element_blank())
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank())
    
    # check performance of real/filled pixels with in situ matchups
    matchups <- dplyr::left_join(df, in_situ_df %>% dplyr::rename(time=time_week), by=c("bin", "time")) %>%
        dplyr::filter(is.finite(ID)) %>%
        dplyr::rename(Real=var, Filled=validation_pixels) %>%
        tidyr::pivot_longer(cols=c(Real, Filled), names_to="Satellite", values_to="sat_chla", values_drop_na=TRUE) %>%
        dplyr::arrange(Satellite, time)
    isp <- matchup_plot(matchups)
    
    return(list(map1=map1,map2=map2,isp=isp,p1=p1,p2=p2,p=p,attr_df=attr_df,week_rmses=week_rmses))
}


# display number of EOFs used to fill the satellite data, and some rmses
display_rmse <- function(attr_df, weeks, week_rmses) {
    cat("Number of EOF:", attr_df$eof, "\n",
        "Total RMSE:", attr_df$rmse, "\n",
        "Week", weeks[1], "RMSE:", week_rmses[1], "\n",
        "Week", weeks[2], "RMSE:", week_rmses[2])
}
