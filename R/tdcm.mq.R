#' @title Utility function to estimate TDCM with multiple Q-matrices.
#'
#' @param data item response data
#' @param qmatrix_list a list of \eqn{I \times A} matrices indicating which
#'  items measure which attributes.
#' @param time.invariance logical argument. If `TRUE` (default), item parameters
#'  are assumed to be equal for all time points. If `FALSE`, item parameters
#'  are not assumed to be equal for all time points.
#' @param rule specific DCM to estimate
#' @param linkfct specific link function for the LCDM
#' @param forget.att vector of attributes not allowing regression/forgetting
#' @param anchor anchor items specified as pairs in a vector
#' @param group A required ``vector`` of integer group identifiers for
#'    multiple group estimation.
#' @param group_invariance logical argument. If `TRUE` (default), item
#'  parameters are assumed to be equal for all groups. If `FALSE`, item
#'  parameters are not assumed to be equal for all groups.
#' @param gdina_extra list of extra arguments to be passed to the gdina function.
#'
#' @keywords internal
#' @noRd
tdcm.mq <- function(
  data,
  qmatrix_list,
  time.invariance = TRUE,
  rule = "LCDM",
  linkfct = "logit",
  forget.att = NULL,
  anchor = NULL,
  group = NULL,
  group_invariance = TRUE,
  gdina_extra = list()
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
  # SETTING THE DESIGN MATRIX
  # ===========================================================================

  design_matrix <- design_matrix_TDCM(
    qmatrix_list,
    time.invariance = time.invariance,
    anchors = anchor,
    rule = rule,
    group = group,
    group_invariance = group_invariance
  )

  # ===========================================================================
  # validate extra args
  # ===========================================================================

  # TODO: Improve error message
  # users can't overwrite the ones we are passing because R throws an error
  # when doing that

  # We get the names for all the arguments from the 'gdina' function
  allowed_args <- setdiff(names(formals(CDM::gdina)), "...")

  # Compared them with the ones passed by the user
  invalid_args <- setdiff(names(gdina_extra), allowed_args)


  # For each invalid argument we emit an error
  for (arg in invalid_args) {
    rlang::arg_match0(arg, allowed_args, arg_nm = "gdina_args")
  }


  # ===========================================================================
  # run model with the new arguments
  # ===========================================================================

  # We write our base parameters as a list
  fn_args <- list(
    data = data,
    q.matrix = qnew,
    linkfct = linkfct,
    method = "ML",
    progress = FALSE,
    delta.designmatrix = design_matrix,
    skillclasses = reduced_space,
    rule = rule,
    group = group,
    invariance = group_invariance
  )

  # we call gdina with our base parameters and the extra ones.
  tdcm <- do.call(CDM::gdina, c(fn_args, gdina_extra))

  return(tdcm)
}
