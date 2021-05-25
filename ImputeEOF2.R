# Temporarily alter metR::ImputeEOF()
# https://github.com/eliocamp/metR/blob/master/R/ImputeEOF.R

# Changes:
#       - return indices of validation points
#       - add temporary simplified/faster code to make a matrix from a datatable
#       - CV_pixels: when you fill gaps using only one year of data, a vector of cross-validation pixels is created.
#           To help compare this to gap filling with multiple years of data, you can re-use that same CV vector for
#           the target year (and more CV pixels will be added for the extra years)

ImputeEOF2 <- function(formula, max.eof = NULL, data = NULL,
                      min.eof = 1, tol = 1e-2, max.iter = 10000,
                      validation = NULL, verbose = interactive(), num_rows,
                      CV_pixels=NULL, all_years, available_time=NULL) {
    
    # check formula and stop if the variable to fill has no missing data
    f <- stringr::str_split(as.character(formula), "~", n = 2)[[1]]
    dcast.formula <- stringr::str_squish(f[stringr::str_detect(f, "\\|")])
    dcast.formula <- as.formula(stringr::str_replace(dcast.formula, "\\|", "~"))
    value.var <- stringr::str_squish(f[!stringr::str_detect(f, "\\|")])
    nas <- sum(is.na(data[[value.var]]))
    if (nas == 0) {
        warning("data has no missing values")
        return(data[[value.var]])
    }
    
    # Build matrix
    
    # This is temporary replacement code specific to DFO gap filling needs (simpler/faster)
    # NOTE1: this won't work if you remove entries with NA "var" from the input dataframe,
    # and possibly not if you are filling multiple years
    all_X <- list()
    all_id <- list()
    for (i in 1:length(all_years)) {
        y <- all_years[i]
        tmp_df <- data %>% dplyr::filter(floor(time)==y)
        all_X[[i]] <- matrix(tmp_df$var, nrow=num_rows)
        tmp_id <- c(matrix(1:nrow(tmp_df), nrow=num_rows))
        if (i>1) {tmp_id <- tmp_id + max(all_id[[i-1]])}
        all_id[[i]] <- tmp_id
    }
    X <- do.call(cbind, all_X)
    id <- as.numeric(do.call(c, all_id))
    
    if (is.null(max.eof)) max.eof <- min(ncol(X), nrow(X))
    gaps <- which(is.na(X))
    
    # get the number of points to use in validation
    if (is.null(validation)) {
        validation <- max(30, 0.1*length(X[-gaps]))
    }
    set.seed(42)    # let's get reproducible in here.
    validation <- sample(seq_along(X)[!seq_along(X) %in% gaps], validation)
    
    
    # use old CV pixels for target year
    if (!is.null(CV_pixels)) {
        inds_bef_target_year <- sum(sapply(1:num_years, function(i) length(available_time[[i]]))) * num_rows
        len_target_year <- length(available_time[[floor(length(available_time)/2)+1]]) * num_rows
        CV_mat <- matrix(CV_pixels, nrow=num_rows)
        CV_pixels <- which(is.finite(CV_mat))
        validation <- validation[validation <= inds_bef_target_year | validation > (inds_bef_target_year+len_target_year)]
        validation <- sort(c(validation, CV_pixels+inds_bef_target_year))
    }
    
    
    eofs <- c(0, min.eof:max.eof)
    X.rec <- X
    # First try, imput with mean or something. Rmse is infinite.
    fill <- mean(X[!is.na(X)])
    X.rec[c(gaps, validation)] <- fill
    rmse <- sqrt(mean((X[validation] - X.rec[validation])^2))
    
    prev <- NULL
    for (i in 2:length(eofs)) {
        # After first guess, impute gaps and validation.
        X.rec <- .ImputeEOF1(X.rec, c(gaps, validation), eofs[i],
                             tol = tol, max.iter = max.iter,
                             verbose = verbose, prev = prev)
        
        prev <- X.rec$prval
        X.rec <- X.rec$X.rec
        
        rmse <- c(rmse, sqrt(mean((X[validation] - X.rec[validation])^2)))
        
        if (verbose == TRUE) {
            cat("\r", "With", eofs[i], "eof - rmse = ", rmse[i])
        }
        
        # Break the loop if we are over the minimum eof asked and, either
        # this rmse is greater than the previous one or current rmse is inf.
        if (rmse[i - 1] - rmse[i] < tol) {
            break
        }
    }
    
    validation_pixels <- matrix(nrow=nrow(X.rec), ncol=ncol(X.rec))
    validation_pixels[validation] <- X.rec[validation]
    
    # Select best eof and make proper imputation.
    eof <- eofs[which.min(rmse)]
    X[gaps] <- fill
    X.rec <- .ImputeEOF1(X, gaps, eof, tol = tol, max.iter = max.iter,
                         verbose = verbose)$X.rec
    
    if (is.data.frame(data)) {
        X.rec <- c(X.rec)[order(id)]
        validation_pixels <- c(validation_pixels)[order(id)]
    }
    
    attr(X.rec, "eof") <- eof
    attr(X.rec, "rmse") <- min(rmse)
    
    return(list(X_filled=X.rec, validation_pixels=validation_pixels))
    
}

.ImputeEOF1 <- function(X, X.na, n.eof, tol = 1e-2, max.iter = 10000,
                        verbose = TRUE, prev = NULL) {
    X.rec <- X
    v <- NULL
    rmse <- Inf
    for (i in 2:max.iter) {
        if (requireNamespace("irlba", quietly = TRUE)) {
            set.seed(42)
            prval <- irlba::irlba(X.rec, nv = n.eof, v = prev)
        } else {
            prval <- base::svd(X.rec, nu = n.eof, nv = n.eof)
            prval$d <- prval$d[1:n.eof]
        }
        v <- prval$v
        R <- prval$u%*%diag(prval$d, nrow = n.eof)%*%t(v)
        rmse <- c(rmse, sqrt(mean((R[X.na] - X.rec[X.na])^2)))
        
        if (rmse[i-1] - rmse[i] > tol) {
            X.rec[X.na] <- R[X.na]
        } else {
            return(list(X.rec = X.rec, prval = prval))
        }
    }
}
