Modified to financial model, have not modified constrained or subj rating


Running:
fit_parameters_mdm: the script to run
fit_parameters_mdm_rating: use subjective rating to fit the SV model 

Functions:
fit_ambigNrisk_model / fit_ambigNrisk_model_Constrained: fitting, separate between unconstrained and constrained models
choice_prob_ambigNrisk: choice probability function
ambig_utility: utility function
getSubjectsInDir: get all the subjects number from the data file folder
exportfig: export figure

After running fit_parameters:
print_fitpar_files: take parameters files and print data to Excel. Including model fitting parameters and choice data.
print_fitpar_files_rating: take parameters files after fitting using subjective rating and print data to Excel. Only model fitting parameters.