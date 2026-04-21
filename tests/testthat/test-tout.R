test_that("printing works (binary)", {
  x <- tout_design(rho_0=0.5, rho_1=0.7, alpha_nom=0.05, beta_nom=0.1, gamma_nom = 0.9)
  expect_snapshot(x)
})
  
test_that("printing works (continuous)", {
  x <- tout_design(rho_0=0, rho_1=0.4, alpha_nom=0.05, beta_nom=0.1, gamma_nom = 0.6, sigma=1)
  expect_snapshot(x)
})


test_that("plotting works (binary)", {
  skip_on_ci()
  x <- tout_design(rho_0=0.5, rho_1=0.7, alpha_nom=0.05, beta_nom=0.1, gamma_nom = 0.9)
  bin_fig <- function() plot(x)
  vdiffr::expect_doppelganger("Binary tout plot", bin_fig)
})

test_that("plotting works (continuous)", {
  skip_on_ci()
  x <- tout_design(rho_0=0, rho_1=0.4, alpha_nom=0.05, beta_nom=0.1, gamma_nom = 0.6, sigma=1)
  cont_fig <- function() plot(x)
  vdiffr::expect_doppelganger("Continuous tout plot", cont_fig)
})

