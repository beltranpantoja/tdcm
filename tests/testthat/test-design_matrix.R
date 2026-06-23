test_that("base design matrix works", {
  qmatrix <- diag(2)
  qmatrix[2, 1] <- 1

  expect_equal(
    create_base_design_matrix(qmatrix, rule = "GDINA"),
    diag(6) # 2 intercepts, 1+2 main effects and 1 interaction
  )
})


test_that("1-PLCDM rule works - 1 attribute", {
  qmatrix <- matrix(1, ncol = 1, nrow = 3)

  # 1 main effect + 3 intercepts
  onepl_design_matrix <- matrix(c(
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 1, 0, 0,
    0, 0, 0, 1,
    0, 1, 0, 0
  ), byrow = TRUE, ncol = 4)

  expect_equal(
    create_base_design_matrix(qmatrix, rule = "1PLCDM"),
    onepl_design_matrix
  )
})

test_that("1-PLCDM rule works - several attributes", {
  qmatrix <- matrix(c(
    1, 0,
    1, 0,
    0, 1
  ), byrow = TRUE, ncol = 2)

  # 2 main effects + 3 intercepts
  onepl_design_matrix <- matrix(c(
    1, 0, 0, 0, 0,
    0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 1, 0, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1
  ), byrow = TRUE, ncol = 5)


  expect_equal(
    create_base_design_matrix(qmatrix, rule = "1PLCDM"),
    onepl_design_matrix
  )
})
