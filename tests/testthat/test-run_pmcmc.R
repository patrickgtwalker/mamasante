library(tidyr)
library(dplyr)
library(readr)
tanz_hist_tx_2015 <- read_csv('Q:/anc_pmcmc/tanz/MAP_TZ_treatment_2020.csv') %>%
  filter(Metric=='Effective Treatment'&Year==2015 & !(Name %in% c('Kaskazini Pemba','Kaskazini Unguja','Kusini Pemba','Kusini Unguja','Mjini Magharibi')))%>%
  mutate(tx_prop=Value/100)%>%
  mutate(Name = ifelse(Name=='Dar-es-salaam','Dar Es Salaam',Name))
tanz_hist_tx_2015_list <- tanz_hist_tx_2015$tx_prop
names(tanz_hist_tx_2015_list) <- tanz_hist_tx_2015$Name
test_run_sifter <- run_pmcmc(data_raw = data_sim, #I've added data_sim to the package for an easy test
                             target_prev = 0.3,
                             n_particles = 20,
                             proposal_matrix = diag(0.5,2),
                             max_param=125,
                             n_steps = 1,
                             n_threads = 1,
                             prop_treated = 0.4,
                             # lag_rates = 10,
                             country = 'Burkina Faso',
                             admin_unit = 'Cascades',
                             seasonality_on = 0,
                             state_check = 0,
                             seasonality_check = 0,
                             start_pf_time = 80,
                             comparison = 'pg',
                             initial = 'informed')
test_that("Testing data processing in run_pmcmc", {

  expect_equal(2 * 2, 4)
})
