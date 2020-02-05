# cancer_cohort_trends
R code for Yang W, Kehm RD, Terry MB. Survival model methods for analyses of cancer incidence trends in young adults. Statistics in Medicine. 2020;1–14. https://doi.org/10.1002/sim.8458

This model code is used for analyses in manuscript Yang W, Kehm RD, and Terry MB 2019 (Survival model methods for analyses of cancer incidence trends in young adults)
CAUTION: for use in other research, please read the code carefully and modify it to suit your own need.

This folder includes the code to generate the simulated data shown in Fig 1 (simulation_to_test_surv.modeling_method.R)
 and all the model code (Functions_to_analyze_cancer_trends_by_cohort.R) to 
	1) fit the age-specific, cohort specific incidence rates to different survival models
 	2) based on the survival model fits, compute the mean hazard rate and use regression to test cohort effect
 	3) make predictions using the combined model from 1 and 2.

To run the code, R needs to be installed first along with the specified packages in the script 
To install the packages altogether, use the script: install.packages

To run the simulations: use the script: simulation_to_test_surv.modeling_method.R

If you would like to use the functions from Functions_to_analyze_cancer_trends_by_cohort.R for your own data, 
 please make sure you organize the data as specified in the model input
