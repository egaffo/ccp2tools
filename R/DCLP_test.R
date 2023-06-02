#' Test circular-to-linear proportion variation
#'
#' Fit and compare Beta binomial models
#'
#' @param id circRNA id to test
#' @param circ matrix with backsplice junction read counts
#' @param lin matrix with the counts of the linearly spliced reads on backsplice
#'   junction
#' @param group data.frame with samples covariates
#' @param full character with the full model formula
#' @param red character with the reduced model formula
#'
#' @return numeric P-value
#' @import aod
#'
#' @examples
test_betabin <- 
  function(id, circ, lin, group, full, red) {
    
    circ <- as.numeric(circ[id, ])
    lin <- as.numeric(lin[id, ])
    tot <- round(lin + circ)
    
    # if there is 0 in the total count vector, the model will fail. So permute 0 to 1
    if ( 0 %in% tot ) {
      tot[tot == 0] <- 1
    }
    
    # Construct data frame
    testdat <- data.frame(tot, circ, group)
    
    ## do test
    # Null model
    red.mod <- as.formula(paste0("cbind(circ, tot - circ) ", red))
    fitNull <- aod::betabin(red.mod, ~ 1, data = testdat)
    
    # Alternative model
    full.mod <- as.formula(paste0("cbind(circ, tot - circ) ", full))
    fitAlt <- aod::betabin(full.mod, ~ 1, data = testdat)
    
    # test models
    a <- aod::anova(fitNull, fitAlt)
    p.value <- a@anova.table[, 11][2]
    p.value
  }

#' Differential circular-to-linear proportion tests
#' 
#' Differential circular-to-linear proportion tests run in parallel
#'
#' @param circ matrix with backsplice junction read counts
#' @param lin matrix with the counts of the linearly spliced reads on backsplice
#'   junction
#' @param groups data.frame with samples covariates
#' @param design character with the model formula
#' @param ... additional parameters passed to BiocParallel::bplapply()
#'
#' @return a named vector of P-values
#' @importFrom BiocParallel bplapply
#'
#' @examples \dontrun{
#' # set BiocParallel parameters for parallel execution
#' biocParams <- BiocParallel::SnowParam(workers = 8, "SOCK", 
#'                                       progressbar = FALSE)
#' 
#' # consider a design formula including batch effects
#' Groups <- data.frame(samples = c("sample1", "sample2", "sample3", "sample4"),
#'                      batch = c("batch1", "batch1", "batch2", "batch2"),
#'                      condition = c("WT", "Treated", "WT", "Treated"))
#' design <- "~ batch + condition"
#' 
#' # get inputs: BJR read counts and linearly spliced reads on the BJ
#' prjs <- c("/home/user/circompara2_batch1",
#'           "/home/user/circompara2_batch2")
#'            
#' ## merge the projects           
#' combs <- combine_ccp2_runs(prjs)
#'
#' ## circRNA BJR matrix 
#' Circ <- as.matrix(data.frame(combs$circ_read_count_mt,
#'                              row.names = "circ_id")) 
#' Linear <- as.matrix(data.frame(combs$lin_read_count_mt, 
#'                                row.names = "circ_id"))
#' 
#' # compute differential CLPs
#' DCLP_pvals <- 
#'   DCLP_test(circ = Circ, 
#'            lin = Linear, 
#'            groups = Groups, 
#'            design = design, 
#'            BPPARAM = biocParams)
#'             
#' # compute adjusted P-values
#' DCLP_padj <- p.adjust(p = DCLP_pvals, method = "BH")
#' }
DCLP_test <- function(circ, lin, groups, design = "~ 1", ...) {
  
  ## compose the reduced formula
  covariates <- strsplit(x = design, "\\+")[[1]]
  if (length(covariates) > 1) {
    red <- paste(covariates[1:length(covariates) - 1], collapse = "+")
  } else {
    red <- "~ 1"  
  }
  
  message("Full model: ", design)
  message("Reduced model: ", red)
  message("Testing ", length(rownames(circ)), " circRNAs ...")
  
  ## fit betabinomial parameters and test in parallel
  ct.pvals <- 
    BiocParallel::bplapply(X = rownames(circ), 
                           FUN = test_betabin,
                           circ = circ, 
                           lin = lin[, colnames(circ)],
                           group = groups[colnames(circ), ], 
                           full = design, 
                           red = red,
                           ...)
  
  setNames(object = unlist(ct.pvals), nm = rownames(circ))
}
