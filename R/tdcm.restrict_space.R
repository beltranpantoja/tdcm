#' Utility function for `restrict_space_matrix()`
#'
#' @keywords internal
#' @noRd
single_restrict_space_matrix <- function(profile_space, col_idx) {
  # We just work with the relevant submatrix for simplicity
  sub_matrix <- profile_space[, col_idx, drop = FALSE]
  last_col <- ncol(sub_matrix)

  # We shift column-wise to compare each stage with the previous one
  comparison_mask <- sub_matrix[, -1] >= sub_matrix[, -last_col]

  # We get the indices of the ones go only go up or stay (i.e. not forgot)
  idx <- rowSums(comparison_mask) == last_col - 1

  # We return a subset of the original profile space passed
  profile_space[idx, ]
}


#' Constraints profiles that can be forgotten
#'
#' Given a full profile space, it constraints some attributes so they can only
#'  increase or stay the same (i.e. students can't forget).
#'
#' @param full_space a matrix that represents all the possible profiles
#' @param forget_attrs which attributes should be restricted (i.e. not forgotten)
#' @param num_attrs how many attributes are there really
#'
#' @keywords internal
#' @noRd
constraints_forget_profiles <- function(full_space, forget_attrs, num_attrs) {
  if (ncol(full_space) %% num_attrs != 0) {
    stop("num_attrs argument is not coherent with the full_space.")
  }

  total_columns <- ncol(full_space)
  restricted_space <- full_space

  for (attr in forget_attrs) {
    # Getting the column indeces
    col_idx <- seq(attr, total_columns, by = num_attrs)
    restricted_space <- single_restrict_space_matrix(restricted_space, col_idx)
  }

  # Return
  restricted_space
}
