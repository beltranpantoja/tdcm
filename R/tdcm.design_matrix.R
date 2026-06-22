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
  # We create base design matrices for all Q-matrices
  design_matrices <- lapply(
    qmatrix_list,
    function(qmat) {
      create_base_design_matrix(qmat, rule)
    }
  )


  # ======================================================
  # INVARIANCE LOGIC
  # if invariant it basically means all items are anchors
  # so the easiest way is to just bind the design matrices
  # (this assumes the qmatrices in the list are all the same)
  # ======================================================

  if (invariance == TRUE) {
    design_matrix <- do.call(rbind, design_matrices)

    # Early return
    return(design_matrix)
  }

  # ======================================================
  # If not invariant, we need to collate them diagonally
  # ======================================================


  # Then we join it into a single design matrix
  design_matrix <- Matrix::.bdiag(design_matrices)


  # If there are no anchors then we just return the design matrix
  if (is.null(anchors)) {
    return(design_matrix)
  }

  # ======================================================
  # ANCHORING LOGIC
  # ======================================================

  # We create the table of the parameters to guide the anchoring process

  qnew <- Matrix::.bdiag(qmatrix_list)
  data <- matrix(
    sample(c(0, 1), 1000 * ncol(qnew), replace = TRUE),
    ncol = ncol(qnew)
  )
  c0 <- TDCM:::tdcm.base(data, qnew, rule)$coef


  # The anchoring matrix established the anchor pairs
  anchor_matrix <- matrix(anchor, ncol = 2, byrow = TRUE)

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
