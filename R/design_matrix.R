#' Create design matrix
#'
#' It creates a design matrix based on the passed rule. This function is useful
#'  for the 1PLCDM implementation and the TDCM function and is not really meant
#'  to be used as a standalone.
#'
#' @param qmatrix a required \eqn{I \times A} matrix indicating which items
#'  measure which attributes.
#' @param rule A string or a vector of itemwise condensation rules that
#'  specific DCM to estimate.
#'
#' @returns A design matrix to be used as `delta.designmatrix` parameter when
#'  calling [CDM::gdina()].
#'
#' @keywords internal
#' @noRd
create_base_design_matrix <- function(qmatrix, rule) {
  if (rule == "1PLCDM") {
    design_matrix <- design_matrix_1plcdm(qmatrix)
  } else {
    # Fake response data to pass to `tdcm.base`
    data <- matrix(
      sample(c(0, 1), 1000 * nrow(qmatrix), replace = TRUE),
      ncol = nrow(qmatrix)
    )

    # This creates a dummy model that we can use to handle the
    # number of parameters that come from the `rule` argument.
    m0 <- TDCM:::tdcm.base(data, qmatrix, rule)

    # Base design matrix
    design_matrix <- diag(nrow(m0$coef))
  }

  # Return
  design_matrix
}


#' Create design matrix for the 1PLCDM model
#'
#' It creates a design matrix for the 1PLCDM model where the main effects are
#'  restricted.
#'
#' @inheritParams create_base_design_matrix qmatrix
#'
#' @inheritSection create_base_design_matrix returns
#'
#' @keywords internal
#' @noRd
design_matrix_1plcdm <- function(qmatrix) {
  # Check properly built
  if (any(rowSums(qmatrix) > 1)) {
    stop("All items should be simple structure.")
  }

  # Make a design matrix for each attribute
  simple_design_matrices <- lapply(
    colSums(qmatrix),
    FUN = function(num_items) {
      # The general design matrix
      design_matrix <- diag(num_items * 2)

      # We get all the indices of the main effects.
      main_effects_idx <- seq(2, nrow(design_matrix), by = 2)

      # We make the main effects be all equal to the first main effect.
      design_matrix[main_effects_idx, ] <- 0
      design_matrix[seq(2, nrow(design_matrix), by = 2), 2] <- 1

      remove_cols <- which(colSums(design_matrix) == 0)

      # Remove empty columns.
      if (length(remove_cols) > 0) {
        design_matrix <- design_matrix[, -remove_cols, drop = FALSE]
      }

      # Return without the empty columns.
      design_matrix
    }
  )

  # Join them diagonally
  as.matrix(
    Matrix::bdiag(simple_design_matrices)
  )
}
