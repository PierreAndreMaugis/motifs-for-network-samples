# motifs-for-network-samples
Code associated with https://arxiv.org/abs/1701.00505

Therein you will find functions (in network_sample_tests.R) to perform the statistical tests described in the companion article. Included is also code to reproduce all the figures (in figures.R), the connectome study (in brain_analysis.R),and the power analyses (in brain_power_experiment.R and bm_power_experiment.R) also presented in the paper.

To perform subgraph counting I use the approach I present in https://arxiv.org/abs/1701.00177 (the code for this is in subgraph_counting_counting.R)

To audit the code: all the graphs being counted are plottable in commented-out code at the top of network_sample_tests.R and I use consistent names across functions; tests are provided in tests.R, which focus on verifying that functions return and plot valid outputs, and that the p-values are uniformly distributed under the null and not so otherwise; the subgraph counting was tested by regressing against igraph’s subgraph counting tools (these test are in the 1701.00177 project; not yet public.)

This code is provided under GNU license. Other license could be provided, upon (unlikely) requests.
Pierre-André Maugis
