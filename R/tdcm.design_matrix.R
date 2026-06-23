#' Create design matrix for TDCM
#'
#' This calls `create_base_design_matrix` and also handles the anchoring.
#'  Invariance is a simple case in which we replicate a single design matrix so
#'  it is handled on the main function.
#'
#' @param qmatrix_list a list of \eqn{I \times A} matrices indicating which
#'  items measure which attributes. This function assumes that if there's a
#'  single qmatrix it was already.
#'  @param invariance logical. If `TRUE` then item parameters will be
#'    constrained to be equal at each time point.
#' @param anchors a list of items to be anchored. They should be given as pairs.
#'  (e.g. `c(1,11,3,15)` would imply that items 1 and 11 are the same and 3 and
#'  3 and 15 too).
#'
#' @inheritParams create_base_design_matrix rule
#'
#' @inheritSection create_base_design_matrix returns
#'
#' @keywords internal
#' @noRd
design_matrix_TDCM <- function(
  qmatrix_list,
  invariance = TRUE,
  anchors = NULL,
  rule = "GDINA"
) {
  # ======================================================
  # COMPLETE QMATRIX
  # ======================================================

  # Normalizing just in case
  qmatrix_list <- lapply(qmatrix_list, as.matrix)
  qnew <- as.matrix(Matrix::.bdiag(qmatrix_list))

  # ======================================================
  # INVARIANCE LOGIC
  # It means we basically have the same design matrix stacked
  # ======================================================

  if (invariance == TRUE) {
    base_design_matrix <- create_base_design_matrix(
      qmatrix_list[[1]],
      rule = rule
    )

    design_matrix_list <- replicate(
      length(qmatrix_list),
      base_design_matrix,
      simplify = FALSE
    )

    design_matrix <- do.call(rbind, design_matrix_list)

    # Early return
    return(design_matrix)
  }

  # ======================================================
  # DESIGN MATRIX
  # ======================================================

  # We create the base design matrix based on the Q-matrix
  design_matrix <- create_base_design_matrix(qnew, rule = rule)


  # If there are no anchors then we just return the design matrix as is
  if (is.null(anchors)) {
    return(design_matrix)
  }

  # ======================================================
  # ANCHORING LOGIC
  # ======================================================

  # We create the table of the parameters to guide the anchoring process
  data <- matrix(
    sample(c(0, 1), 1000 * ncol(qnew), replace = TRUE),
    ncol = ncol(qnew)
  )
  c0 <- TDCM:::tdcm.base(data, qnew, rule)$coef


  # The anchoring matrix established the anchor pairs
  anchor_matrix <- matrix(anchors, ncol = 2, byrow = TRUE)

  # We iterate over each each anchor-pair
  for (row in seq_len(nrow(anchor_matrix))) {
    # items have multiple parameters, so we get the indices for the source and target
    source_item_idx <- which(c0[, "itemno"] == anchor_matrix[row, 1])
    target_item_idx <- which(c0[, "itemno"] == anchor_matrix[row, 2])

    design_matrix[target_item_idx, ] <- design_matrix[source_item_idx, ]
  }


  # The constrained items will correspond to columns with only 0s
  remove_cols <- which(colSums(design_matrix) == 0)


  # Return
  design_matrix[, -remove_cols]
}
