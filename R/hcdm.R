#' Estimation of the Hierarchical Diagnostic Classification Model (HDCM)
#'
#' `hcdm()` estimates the hierarchical diagnostic classification model
#' (HCDM; Templin & Bradshaw, 2025).
#'
#' @inheritParams tdcm data q.matrix linkfct
#'
#' @param hierarchy A string containing restrictions of the hierarchy.
#'  For example, `'Attr1 > Attr2 > Attr3'` means that attribute 3 requires
#'  attribute 2, which requires attribute 1.
#'
#' @param ... Extra arguments to pass to the [CDM::gdina()] call.
#'
#' @references
#' Templin J, Bradshaw L. (2014). Hierarchical Diagnostic Classification Models:
#'  A Family of Models for Estimating and Testing Attribute Hierarchies.
#'  _Psychometrika_, 79_(2), 317-339. <doi:10.1007/s11336-013-9362-0>.
#'
#' @return An object of class `hcdm` which inherits from [CDM::gdina]
#'
#' @export
hcdm <- function(
  data,
  q.matrix,
  hierarchy,
  linkfct = "logit",
  ...
) {
  # We get the reduced skill space
  skillspace <- CDM::skillspace.hierarchy(
    hierarchy,
    skill.names = colnames(q.matrix)
  )

  zeroprob.skillclasses <- skillspace$zeroprob.skillclasses

  # Then we get the new Mj to pass to the gdina model
  # This restricts the parameters to only the ones that are possible
  restricted_Mj <- get_HDCM_Mj(hierarchy, q.matrix)


  # Estimate model
  # The CDM documentation recommends using avoid.zeroprobs=TRUE
  # to avoid algorithmic instabilities
  model <- CDM::gdina(
    data,
    q.matrix,
    zeroprob.skillclasses = zeroprob.skillclasses,
    Mj = restricted_Mj,
    avoid.zeroprobs = TRUE,
    link = linkfct,
    rule = "GDINA", # The rules are enforced through the hierarchy only
    ...
  )

  # =======================================================================
  # CLASS DEFINITION
  # We extend the gdina class so we can define it's own methods and attributes
  # =======================================================================
  class(model) <- c("hcdm", class(model))


  model$hierarchy <- hierarchy
  # =======================================================================

  return(model)
}


#' Builds the design matrix M following a hierarchy
#'
#' @keywords internal
#' @noRd
get_HDCM_Mj <- function(hierarchy_str, q.matrix) {
  # 1. Get the skillspace of the hierarchy
  skillspace <- CDM::skillspace.hierarchy(
    hierarchy_str,
    skill.names = colnames(q.matrix)
  )


  # For each item we get the submatrix of possible profiles
  # we will use this to filter the base design matrix Aj.
  reduced_profiles_item <- lapply(
    seq_len(nrow(q.matrix)),
    FUN = \(i) {
      skillspace$skillspace.reduced[, q.matrix[i, ] == 1, drop = FALSE]
    }
  )

  # Generate fake data to get a model to use as base
  data <- matrix(
    sample(c(0, 1), nrow(q.matrix) * 1000, TRUE),
    ncol = nrow(q.matrix)
  )

  colnames(data) <- paste0("V", seq_len(nrow(q.matrix)))

  mod <- CDM::gdina(data, q.matrix = q.matrix, maxit = 1, progress = FALSE)


  # We check the rows of the A design matrix which has all possible profiles
  # against the reduced profile space
  parameters_idx_by_item <- mapply(
    FUN = \(A, B) {
      A_rows <- apply(A, 1, paste, collapse = "\r")
      B_rows <- apply(B, 1, paste, collapse = "\r")

      A_rows %in% B_rows
    },
    A = mod$Aj,
    B = reduced_profiles_item,
    SIMPLIFY = FALSE
  )

  # Now we have the columns which are related to permisible parameters
  # So we use that to generate a new reduced Mj

  new_Mj <- mapply(
    FUN = \(M, idx) {
      M[[2]] <- M[[2]][idx]
      M[[1]] <- M[[1]][, idx]
      M
    },
    M = mod$Mj,
    idx = parameters_idx_by_item,
    SIMPLIFY = FALSE
  )

  return(new_Mj)
}
