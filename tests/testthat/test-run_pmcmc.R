test_run_sifter <- run_pmcmc(data_raw = data_sim, #I've added data_sim to the package for an easy test
                             init_EIR = 10,
                             n_particles = 20,
                             proposal_matrix = diag(0.5,2),
                             max_param=125,
                             n_steps = 1,
                             n_threads = 1,
                             # lag_rates = 10,
                             country = 'Burkina Faso',
                             admin_unit = 'Cascades',
                             seasonality_on = 0,
                             state_check = 0,
                             seasonality_check = 0,
                             start_pf_time = 80,
                             comparison = 'u5',
                             stoch_param = 'betaa',
                             initial = 'informed')
test_that("Testing data processing in run_pmcmc", {

  expect_equal(2 * 2, 4)
})
