remove_batch_effect_binomial <-
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
                  FUN = function(id, x, mod){
                    stats::glm.fit(x = mod, 
                                   y = as.numeric(x[id, ]),
                                   family = stats::binomial())
                  },
                  x = x,
                  mod = full.mod)
    
    beta <- sapply(fit, 
                   function(f)f$coefficients[-(seq_len(ncol(design)))],
                   simplify = T)
    
    beta[is.na(beta)] <- 0
    
    x <- as.matrix(x)
    LM <- car::logit(x, percents = F, adjust = 0)
    boot::inv.logit(LM - beta %*% t(X.batch))
    
  }
