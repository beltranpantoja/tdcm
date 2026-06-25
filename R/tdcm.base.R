#' Utility function for base estimation in TDCMs.
#'
#' This function allows to leverage the CDM package to handle the rules and
#' other parameters. These can then be used to tweak the design matrices and
#' similar objects.
#'
#' @param data item response data
#' @param q.matrix.induced The induced Q-matrix.
#' @param rule A string or a vector of itemwise condensation rules that specific DCM to estimate.
#'
#' @return An object of class `gdina` returned by the internal call to [CDM::gdina()].
#'
#' @keywords internal
#' @noRd
tdcm.base <- function(
  data,
  q.matrix.induced,
  rule,
  skillclasses = NULL,
  group = NULL,
  group_invariance = TRUE
) {
  tdcm.1 <- suppressWarnings(CDM::gdina(
    data,
    q.matrix.induced,
    rule = rule,
    linkfct = "logit", # use logit link function for all items
    method = "ML", # directly maximize the log-likelihood function
    mono.constr = TRUE, # meet monotonicity constraints in estimation
    progress = FALSE, # do NOT print the iteration progress
    maxit = 1, # We just want the base constructions
    skillclasses = skillclasses,
    group = group,
    invariance = group_invariance
  )) # tdcm.1
  return(tdcm.1)
} # tdcm.base
