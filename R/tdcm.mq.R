#' @title Utility function to estimate TDCM with multiple Q-matrices.
#'
#' @param data item response data
#' @param qmatrix_list a list of \eqn{I \times A} matrices indicating which
#'  items measure which attributes.
#' @param invariance invariance assumption (T or F)
#' @param rule specific DCM to estimate
#' @param linkfct specific link function for the LCDM
#' @param forget.att vector of attributes not allowing regression/forgetting
#' @param anchor anchor items specified as pairs in a vector
#' @keywords internal
#' @noRd
tdcm.mq <- function(
  data,
  qmatrix_list,
  invariance = TRUE,
  rule = "LCDM",
  linkfct = "logit",
  forget.att = NULL,
  anchor = NULL
) {
  # ===========================================================================
  # SETTING THE COMPLETE QMATRIX
  # ===========================================================================

  # We convert the elements to matrices
  qmatrix_list <- lapply(
    qmatrix_list,
    as.matrix
  )

  # We join diagonally and convert to a regular R matrix
  qnew <- as.matrix(Matrix::.bdiag(qmatrix_list))

  # ===========================================================================
  # SETTING THE DESIGN MATRIX
  # ===========================================================================

  design_matrix <- design_matrix_TDCM(
    qmatrix_list,
    invariance = invariance,
    anchors = anchor,
    rule = rule
  )

  # ===========================================================================
  # (NOT) FORGET ATTRIBUTE LOGIC
  # ===========================================================================


  # This allows us to get the complete profile space
  m0 <- tdcm.base(
    data,
    q.matrix.induced = qnew,
    rule = rule
  )

  # Default is no restriction
  reduced_space <- m0$attribute.patt.splitted

  # We constraint the space
  if (!is.null(forget.att)) {
    reduced_space <- constraints_forget_profiles(
      full_space = reduced_space,
      forget_attrs = forget.att,
      num_attrs = ncol(qmatrix_list[[1]])
    )
  }

  # ===========================================================================
  # run model with the new arguments
  # ===========================================================================

  tdcm <- CDM::gdina(
    data,
    qnew,
    linkfct = linkfct,
    method = "ML",
    progress = FALSE,
    delta.designmatrix = design_matrix,
    skillclasses = reduced_space,
    rule = rule
  )

  return(tdcm)
}
