# msclust: monotonic statistical clustering

Code to supplement the manuscript entitled "Data-driven clustering of infectious disease incidence into age groups". 

Authors: Rami Yaari, Amit Huppert and Itai Dattner.

## Content

- util.R: various utility functions.

- run_sir.R : generates simulated age group incidence data using SIR model .

- msclust.R: implementation of the msclust algorithm.

- msclust_test.R: testing the peformance of the msclust algorithm using Monte-Carlo simulations of age-group incidence data.

- msclust_test_with_bagging.R: testing the peformance of the msclust+ algorithm (msclust with bagging) using Monte-Carlo simulations of age-group incidence data.

- gmm_test.R: testing the peformance of the Gaussian mixture model algorithm using Monte-Carlo simulations of age-group incidence data.

- sigclust2_test.R: testing the peformance of the sigclust2 algorithm using Monte-Carlo simulations of age-group incidence data.

- run_alpha_tests.R: running tests to assess the type-I error rate of the various algorithms.

- run_power_tests.R: running tests to assess the power of the various algorithms.

- covid19_age_group_clustering.R: runnning the msclust+ algorithm on real incidence data of Covid-19 from Israel, to obtain age-group clustering per epidemic wave and sector of the population. The data used by this code cannot be made available at the current time for privacy considerations. 

#### Support and Contributions:
For support and bug reports send an email to: ramiyaari@gmail.com or open an issue [<a href="https://github.com/ramiyaari/simode/issues">here</a>].
