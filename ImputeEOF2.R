# Temporarily alter metR::ImputeEOF()
# https://github.com/eliocamp/metR/blob/master/R/ImputeEOF.R

# Changes:
#       - return indices of validation points
#       - add temporary simplified/faster code to make a matrix from a datatable

ImputeEOF2 <- function(formula, max.eof = NULL, data = NULL,
                      min.eof = 1, tol = 1e-2, max.iter = 10000,
                      validation = NULL, verbose = interactive(), num_rows) {
    
    # Build matrix.
    f <- stringr::str_split(as.character(formula), "~", n = 2)[[1]]
    dcast.formula <- stringr::str_squish(f[stringr::str_detect(f, "\\|")])
    dcast.formula <- as.formula(stringr::str_replace(dcast.formula, "\\|", "~"))
    value.var <- stringr::str_squish(f[!stringr::str_detect(f, "\\|")])
    nas <- sum(is.na(data[[value.var]]))
    if (nas == 0) {
        warning("data has no missing values")
        return(data[[value.var]])
    }
    
    # # temporary replacement code specific to dfo gap filling needs (simpler/faster)
    # # note: this won't work if you remove entries with NA "var" from the input dataframe,
    # # and possibly not if you are filling multiple years
    # # ALSO NOTE: I don't understand why, but byrow=TRUE is the correct way to do this here...
    # #     even though df should be written column-wise to get row_pixel x column_day
    # id = c(matrix(1:nrow(data), nrow=num_rows, byrow=TRUE))
    # X = matrix(data$var, nrow=num_rows, byrow=TRUE)
    
    g <- .tidy2matrix(data, dcast.formula, value.var)
    data$ff19bdd67ff5f59cdce2824074707d20 <- 1:nrow(data)
    id <- c(.tidy2matrix(data, dcast.formula, "ff19bdd67ff5f59cdce2824074707d20")$matrix)
    X <- g$matrix
    
    if (is.null(max.eof)) max.eof <- min(ncol(X), nrow(X))
    gaps <- which(is.na(X))
    # if (length(gaps) == 0) return(X)

    # get the number of points to use in validation
    if (is.null(validation)) {
        validation <- max(30, 0.1*length(X[-gaps]))
    }
    set.seed(42)    # let's get reproducible in here.
    validation <- sample(seq_along(X)[!seq_along(X) %in% gaps], validation)

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
    validation_pixels <- c(validation_pixels)[order(id)]
    
    # Select best eof and make proper imputation.
    eof <- eofs[which.min(rmse)]
    X[gaps] <- fill
    X.rec <- .ImputeEOF1(X, gaps, eof, tol = tol, max.iter = max.iter,
                         verbose = verbose)$X.rec
    
    if (is.data.frame(data)) {
        X.rec <- c(X.rec)[order(id)]
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

.tidy2matrix <- function(data, formula, value.var, fill = NULL, ...) {
    row.vars <- all.vars(formula[[2]])
    col.vars <- all.vars(formula[[3]])
    data <- data.table::as.data.table(data)
    data[, row__ := .GRP, by = c(row.vars)]
    data[, col__ := .GRP, by = c(col.vars)]
    if (is.null(fill)){
        fill <- 0
        # rowdims <- data[col__ == 1, (row.vars), with = FALSE]
        # coldims <- data[row__ == 1, (col.vars), with = FALSE]
    } else {
        # rowdims <- unique(data[, (row.vars), with = FALSE])
        # coldims <- unique(data[, (col.vars), with = FALSE])
    }
    rowdims <- unique(data[, (row.vars), with = FALSE])
    coldims <- unique(data[, (col.vars), with = FALSE])
    data.m <- matrix(fill[1], nrow = max(data[["row__"]]),
                     ncol = max(data[["col__"]]))
    data.m[cbind(data[["row__"]], data[["col__"]])] <- data[[value.var]]
    
    return(list(matrix = data.m,
                coldims = coldims,
                rowdims = rowdims))
}
