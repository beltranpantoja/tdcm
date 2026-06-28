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
  rule = "LCDM",
  linkfct = "logit",
  time.invariance = TRUE,
  forget.att = NULL,
  anchor = NULL,
  group = NULL,
  group_invariance = TRUE,
  hierarchy = NULL,
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

  # This names are just to avoid errors. Then the summary should use the
  # original ones
  colnames(qnew) <- paste0("Att", seq_len(ncol(qnew)))

  # ===========================================================================
  # HANDLING HIERARCHY OF ATTRIBUTES AND PROFILE SPACE
  # ===========================================================================

  # By default the M matrix is NULL and the profile space is complete
  Mj <- NULL
  profile_space <- CDM::skillspace.full(colnames(qmatrix_list[[1]]))

  # If there's a hierarchy then we redefine these values
  if (!is.null(hierarchy)) {
    Mj <- get_HDCM_Mj(
      hierarchy,
      q.matrix = do.call(rbind, qmatrix_list)
    )

    profile_space <- CDM::skillspace.hierarchy(
      hierarchy,
      skill.names = colnames(qmatrix_list[[1]])
    )$skillspace.reduced
  }


  # Now we extend the profile space to fill all the time points
  time_points <- length(qmatrix_list)
  profile_space <- extend_profile_space(profile_space, time_points)

  # ===========================================================================
  # (NOT) FORGET ATTRIBUTE LOGIC
  # We limit the rows depending on which attributes can't be forgotten.
  # ===========================================================================

  if (!is.null(forget.att)) {
    profile_space <- constraints_forget_profiles(
      full_space = profile_space,
      forget_attrs = forget.att,
      num_attrs = ncol(qmatrix_list[[1]])
    )
  }

  # ===========================================================================
  # DETERMINING THE ZERO PROB PROFILES
  # we need to pass the vector of not possible profiles following the CDM order
  # ===========================================================================

  # We generate all the possible profiles, no restrictions
  full_profile_space <- CDM::skillspace.full(colnames(qnew))


  # This is code just iterates over the complete space and see which ones
  # are present in our reduced space
  full_space_rows <- apply(full_profile_space, 1, paste, collapse = "")
  reduced_space_rows <- apply(profile_space, 1, paste, collapse = "")

  # We want to know which ones are NOT present. Those are the zeroprob classes
  zeroprob_idx <- which((full_space_rows %in% reduced_space_rows) == FALSE)

  # ===========================================================================
  # SETTING THE DESIGN MATRIX
  # ===========================================================================

  design_matrix <- design_matrix_TDCM(
    qmatrix_list,
    time.invariance = time.invariance,
    anchors = anchor,
    rule = rule,
    group = group,
    group_invariance = group_invariance,
    Mj = Mj
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
    rule = rule,
    group = group,
    Mj = Mj,
    zeroprob.skillclasses = zeroprob_idx,
    avoid.zeroprobs = TRUE,
    reduced.skillspace = FALSE,
    invariance = group_invariance
  )

  # we call gdina with our base parameters and the extra ones.
  tdcm <- do.call(CDM::gdina, c(fn_args, gdina_extra))

  return(tdcm)
}


#' Utility function to get the complete profile space over time.
#'
#' It works by chaining all possible combinations of the profile space n times.
#'
#' @keywords internal
#' @noRd
extend_profile_space <- function(profile_space, n) {
  row_indices <- replicate(
    n,
    seq_len(nrow(profile_space)),
    simplify = FALSE
  )

  idx_combinations <- expand.grid(row_indices)

  result <- do.call(cbind, lapply(seq_len(n), function(i) {
    profile_space[idx_combinations[[i]], , drop = FALSE]
  }))
}
