#' Remove Batch Effects by fitting a Qquasi-binomial GLM
#'
#' A function to remove batch effects from a matrix of proportion values. Like
#' the homonymous function from the limma package but for values between 0 and 1.
#' Useful to correct batch effects from CLP matrices.
#'
#' @author Enrico Gaffo
#' @param x numeric matrix containing proportion (i.e., real numbers between 0
#'  and 1) values for a series of samples. Rows correspond to probes and
#'  columns to samples.
#' @param batch factor or vector indicating batches.
#' @param covariates matrix or vector of numeric covariates to be adjusted for.
#' @param design model design matrix
#' @param ... not used
#'
#' @import stats car boot
#'
#' @export
#'
#' @seealso limma::removeBatchEffect()
#' @return A numeric matrix of proportion values with batch and covariate
#'   effects removed.
remove_batch_effect_quasibin <-
  function(x,
            batch,
            covariates = NULL,
            design = matrix(1, ncol(x), 1),
            ...) {

    if (is.null(batch))
      return(as.matrix(x))

    if (!is.null(batch)) {
      Batch <- as.factor(batch)
      stats::contrasts(Batch) <- stats::contr.sum(levels(Batch))
      Batch <- stats::model.matrix(~Batch)[, -1, drop = FALSE]
    }

    if (!is.null(covariates))
      covariates <- as.matrix(covariates)

    X.batch <- cbind(Batch, covariates)

    full.mod <- cbind(design, X.batch)

    fit <- lapply(X = rownames(x),
                  FUN = fit_quasibinomial,
                  x = x,
                  mod = full.mod)

    beta <- sapply(fit, function(f)f$coefficients[-(seq_len(ncol(design)))],
                   simplify = T)

    beta[is.na(beta)] <- 0

    x <- as.matrix(x)
    LM <- car::logit(x, percents = F, adjust = 0)
    boot::inv.logit(LM - beta %*% t(X.batch))
  }
