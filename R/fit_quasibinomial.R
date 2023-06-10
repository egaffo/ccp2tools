#' Fit a GLM of the quasibinomial family
#'
#' Auxiliary function to fit a Generalized Linear Model of the quasibinomial
#' family to a row of a matrix
#'
#' @author Enrico Gaffo
#' @param id an identifier to select the matrix row
#' @param x a numeric matrix
#' @param mod a matrix representing the model to fit. For instance, as returned
#'   by the model.matrix() function
#' @return an object of class inheriting from "glm", as the output of the
#'   glm.fit() function
#'   
#' @import stats
#' @export
#' 
fit_quasibinomial <-
  function(id, x, mod) {

    abund <- as.numeric(x[id, ])

    stats::glm.fit(x = mod, y = abund, family = stats::quasibinomial())
  }
