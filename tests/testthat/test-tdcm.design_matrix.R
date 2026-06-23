test_that("design matrix invariance works", {
  qmatrix <- diag(2)
  qmatrix_list <- list(qmatrix, qmatrix)

  expect_equal(
    design_matrix_TDCM(qmatrix_list, invariance = FALSE),
    diag(8)
  )

  expect_equal(
    design_matrix_TDCM(qmatrix_list, invariance = TRUE),
    rbind(diag(4), diag(4))
  )
})


test_that("design matrix with anchors works", {
  qmatrix <- diag(2)
  qmatrix_list <- list(qmatrix, qmatrix)

  base_design_matrix <- design_matrix_TDCM(qmatrix_list, invariance = FALSE)

  # Anchoring the first and last items will make these two rows the same
  base_design_matrix[c(7, 8), ] <- base_design_matrix[c(1, 2), ]

  # The last two columns dissapear
  base_design_matrix <- base_design_matrix[, -c(7, 8)]

  expect_equal(
    design_matrix_TDCM(qmatrix_list, invariance = FALSE, anchors = c(1, 4)),
    base_design_matrix
  )
})
