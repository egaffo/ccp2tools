#' Remove Batch Effects by fitting a Beta Binomial GLM
#'
#' Similar to the homonymous function from the limma package but the data are
#' fitted with a Beta-binomial distribution.
#'
#' @author Enrico Gaffo
#'
#' @param x a numeric matrix with the BJR counts
#' @param y a numeric matrix with the counts of linearly spliced reads on the
#'   backsplice junctions
#' @param batch a vector indicating the names of the batch variables which
#'   effect has to be removed
#' @param design a vector indicating the names of the variables which effect has
#'   to be kept
#' @param groups a data frame with the sample groups for the batch and main
#'   effect vriables
#' @param ... Further arguments passed to optim within aod::betabin
#'
#' @import aod boot car data.table stats
#' @export
#' 
#' @seealso limma::removeBatchEffect()
#' @return A numeric matrix of proportion values with batch effects removed.
#' 
#' @examples \dontrun{
#' ### test code ####
#' 
#' source("R/combine_ccp2_runs.R")
#' combs <- combine_ccp2_runs(files = c("/sharedfs01/enrico/CLL/analysis/CCP2/",
#'                                      "/sharedfs01/enrico/CLL/analysis/PRJNA432966/"),
#'                            merge_circs = T, merge_lin_bks = T, recycle_existing_lincount = T,
#'                            merge_lin = F)
#' 
#' circs <- as.matrix(data.frame(combs$circ_read_count_mt, row.names = "circ_id", check.names = F))
#' 
#' lins <- as.matrix(data.frame(combs$lin_read_count_mt, row.names = "circ_id", check.names = F))
#' lins <- lins[rownames(circs), colnames(circs)]
#' lins[is.na(lins)] <- 0
#' 
#' colnames(circs) <- sub("_.*", "", colnames(circs))
#' colnames(lins) <- sub("_.*", "", colnames(lins))
#' 
#' batch1_data <- 
#'   data.table::fread("/sharedfs01/enrico/CLL/analysis/R_CLL/data/CLL_meta_short.csv")
#' batch2_data <-
#'   data.table::fread("/sharedfs01/enrico/CLL/analysis/R_CLL/data/PRJNA432966_meta.csv")
#' meta <- data.table::rbindlist(l = list(batch1 = batch1_data,
#'                                        batch2 = batch2_data),
#'                               use.names = T, fill = TRUE, idcol = "Batch")
#' meta <- meta[`SAMPLE ID` != "1224666"]
#' meta[Batch == "batch2", `:=`(CONDITION = "B-CELL", TRANSLOCATION = "NONE")]
#' meta_df <- data.frame(meta, row.names = "SAMPLE ID", check.names = F)
#' meta_df$COND_TRANS <- paste(meta_df$CONDITION, meta_df$TRANSLOCATION, sep = "_")
#'
#' batch <- "Batch"
#' design <- "~ COND_TRANS"
#' 
#' x <- circs[rowSums(circs[, rownames(meta_df)] > 5) >= 9, rownames(meta_df)]
#' y <- lins[rownames(x), colnames(x)]
#' 
#' prps <- x / (x + y)
#' 
#' prps[is.na(prps)] <- 0
#' prps <- prps[rowSums(prps > 0) >= 9, ]
#' prps <- prps[order(matrixStats::rowVars(prps), decreasing = T)[1:100], ]
#' 
#' # PCA originaldata ####
#' pcs <- prcomp(t(prps), center = T)
#' df <- cbind(meta_df, pcs$x[rownames(meta_df), ])
#' 
#' library(ggplot2)
#' ggplot(df, aes(x = PC1, y = PC2, color = COND_TRANS, shape = Batch)) +
#'  geom_point(size = 3) +
#'   ggtitle("original")
#' 
#' ## MDS
#' # ggplot(cbind(cmdscale(dist(t(prps)))[rownames(meta_df), ], meta_df), 
#' #        aes(x = `1`, y = `2`, color = COND_TRANS, shape = Batch)) +
#' #   geom_point(size = 3)
#' 
#' # PCA logit scale data ####
#' pcs <- prcomp(t(car::logit(prps, percents = F, adjust = 1e-6)), 
#'               center = T, scale. = T)
#' df <- cbind(meta_df, pcs$x[rownames(meta_df), ])
#' 
#' ggplot(df, aes(x = PC1, y = PC2, color = COND_TRANS, shape = Batch)) +
#'   geom_point(size = 3) +
#'   ggtitle("original - logit scale")
#' 
#' # remove batch betabin ####
#' rbe <- remove_batch_effect_betabin(x = x[rownames(prps), ], 
#'                                    y = y[rownames(prps), ], 
#'                                    groups = meta_df,
#'                                    batch = "Batch", 
#'                                    design = "~COND_TRANS")
#' # summary(rbe)
#' rbe[is.na(rbe)] <- 0
#' 
#' pcs <- prcomp(t(rbe), center = T)
#' df <- cbind(meta_df, pcs$x[rownames(meta_df), ])
#' 
#' ggplot(df, aes(x = PC1, y = PC2, color = COND_TRANS, shape = Batch)) +
#'   geom_point(size = 3) +
#'   ggtitle("rbe betabin")
#' 
#' # remove batch quasibin ####
#' source("R/remove_batch_effect_quasibin.R")
#' source("R/fit_quasibinomial.R")
#' rbe.qsb <- remove_batch_effect_quasibin(x = prps,
#'                                         batch = meta_df$Batch,
#'                                         design = model.matrix(~COND_TRANS,
#'                                                               data = meta_df))
#' rbe.qsb[is.na(rbe.qsb)] <- 0
#' 
#' pcs <- prcomp(t(rbe.qsb), center = T)
#' df <- cbind(meta_df, pcs$x[rownames(meta_df), ])
#' 
#' ggplot(df, aes(x = PC1, y = PC2, color = COND_TRANS, shape = Batch)) +
#'   geom_point(size = 3) +
#'   ggtitle("rbe quasibin")
#' 
#' # remove batch binomial ####
#' source("R/remove_batch_effect_binomial.R")
#' rbe.bin <- remove_batch_effect_binomial(x = prps,
#'                                         batch = meta_df$Batch,
#'                                         design = model.matrix(~COND_TRANS,
#'                                                               data = meta_df))
#' rbe.bin[is.na(rbe.bin)] <- 0
#' 
#' pcs <- prcomp(t(rbe.bin), center = T)
#' df <- cbind(meta_df, pcs$x[rownames(meta_df), ])
#' 
#' ggplot(df, aes(x = PC1, y = PC2, color = COND_TRANS, shape = Batch)) +
#'   geom_point(size = 3) +
#'   ggtitle("rbe bin")
#' 
#' # remove batch limma ####
#' rbe.lim <- limma::removeBatchEffect(x = prps,
#'                                     batch = meta_df$Batch,
#'                                     design = model.matrix(~COND_TRANS,
#'                                                           data = meta_df))
#' # rbe.lim[is.na(rbe.lim)] <- 0
#' 
#' pcs <- prcomp(t(rbe.lim), center = T, scale. = T)
#' df <- cbind(meta_df, pcs$x[rownames(meta_df), ])
#' 
#' ggplot(df, aes(x = PC1, y = PC2, color = COND_TRANS, shape = Batch)) +
#'   geom_point(size = 3) +
#'   ggtitle("rbe limma")
#' 
#' 
#' # remove batch limma logit ####
#' rbe.limlg <- limma::removeBatchEffect(x = car::logit(prps, percents = F,
#'                                                      adjust = 1e-6),
#'                                       # rbe.limlg <- limma::removeBatchEffect(x = log(prps+1e-6),
#'                                       batch = meta_df$Batch,
#'                                       design = model.matrix(~COND_TRANS,
#'                                                             data = meta_df))
#' # rbe.limlg[is.na(rbe.limlg)] <- 1e-6
#' # rbe.limlg <- boot::inv.logit(rbe.limlg) + 1e-6
#' 
#' pcs <- prcomp(t(rbe.limlg), center = T, scale. = T)
#' # summary(pcs)
#' df <- cbind(meta_df, pcs$x[rownames(meta_df), ])
#' 
#' ggplot(df, aes(x = PC1, y = PC2, color = COND_TRANS, shape = Batch)) +
#'   geom_point(size = 3) +
#'   ggtitle("rbe limma on logit scaled")
#' 
#' 
#' }
remove_batch_effect_betabin <-
  function(x, y,
           groups,
           batch = NULL,
           design = "~ 1",
           ...) {
    
    data <- 
      merge(data.table::melt(data.table::as.data.table(x,
                                                       keep.rownames = "ids"),
                             id.vars = "ids"),
            data.table::melt(data.table::as.data.table(y,
                                                       keep.rownames = "ids"),
                             id.vars = "ids"),
            by = c("ids", "variable"))[, .(ids, variable,
                                           x = value.x, 
                                           n = value.x + value.y)]
    
    data <- data.table::as.data.table(cbind(groups[data$variable, ], data))
    # if there is 0 in the total count vector, the model will fail
    data[n == 0, n := 1]
    
    props <- x / (x + y[rownames(x), colnames(x)])
    props[is.na(props)] <- 0
    
    if (is.null(batch) | gsub(" ", "", design) == "~1") {
      return(props)
    }
    
    # if (!is.null(batch)) {
    #   
    #   X.batch <- 
    #     as.matrix(data.frame(sapply(batch, 
    #                                 function(b, g){
    #                                   X.batch <- as.factor(groups[[b]])
    #                                   stats::contrasts(X.batch) <- stats::contr.sum(levels(X.batch))
    #                                   X.batch <- stats::model.matrix(~X.batch)[, -1, drop = FALSE]
    #                                   colnames(X.batch) <- b
    #                                   X.batch
    #                                 }, 
    #                                 g = groups, 
    #                                 simplify = F, 
    #                                 USE.NAMES = T)))
    #   
    # }
    if (!is.null(batch)) {
      Batch <- as.factor(groups[, batch])
      stats::contrasts(Batch) <- stats::contr.sum(levels(Batch))
      Batch <- stats::model.matrix(~Batch)[, -1, drop = FALSE]
    }
    
    X.batch <- cbind(Batch, NULL)
    
    full_design <- paste(sub("~", "", design), 
                         paste0(batch, collapse = "+"),
                         sep = "+")
    full_formula <- as.formula(paste0("cbind(x, n - x) ~ ", full_design))
    
    fit <- lapply(X = split(x = data, f = data$ids), 
                  FUN = function(d, formula) {
                    aod::betabin(formula = formula, random = ~1,
                                 data = d, ...)
                  },
                  formula = full_formula)
    
    design_modmat <- as.matrix(stats::model.matrix(object = as.formula(design),
                                                   data = groups))
    betas <- 
      sapply(X = fit, 
             FUN = function(f)aod::coef(f)[-(seq_len(ncol(design_modmat)))],
             simplify = T)
    
    betas[is.na(betas)] <- 0
    
    LM <- car::logit(props[names(fit), ], percents = F, adjust = 0)
    # LM <- car::logit(props[names(fit), ], percents = F, adjust = 1e-9)
    # LM <- boot::logit(props)
    boot::inv.logit(LM - betas %*% t(-X.batch))
  }
