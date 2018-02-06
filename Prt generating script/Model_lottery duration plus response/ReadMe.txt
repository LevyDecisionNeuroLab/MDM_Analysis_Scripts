Functions:
PTB_Protocol_Gen_ver(no.)
PTB_Protocol_OnsetExtract
getBolckOrder
getBlockOnsetOffset
getSubjectsInDir
ambig_utility

Run:
PTB_Protocol_Setup_ver(no.)
prt2sdm
add_zeros_columns, add_zeros_columns_ver(no.), adjust_columns_ver(no.)
copysdm2folder

-----------------------------------------------------------------------------------

Versions:
Obosolete ones:


Ver3: 
Got from Lital.
Display+response, each has parametric modulator, 8 or 16 predictors


Ver4:
Outcome magnitude simple binary predictors, seperate for ambiguity and
risk.
- For each domain, 8 seperate binary predictors (4magnitude x
2risk/amb).
- Display and response. So 16 predictors for each domain. 32 in total.


Ver5:
Ver5 create 4 binary predictors for outcome magnitude, do not seperate
for risk/ambiguity, but seperate for each domain.

Ver7:
All others same as ver3, except that response predictor does not have parametric modulator
Ver3 should be obsolete ever after.
8 or 12 predictors
Need to use "add_zeros_columns_ver7"
Include the CV (chosen subjective value) in parametric modulator

Ver6:
Deleted, never used

--------------------------------------------------------------------------------------

In use:

Ver5_1:
Improve from Ver5
4 binary predictors for outcome magnitude, do not seperate for risk/ambiguity, but seperate for each domain.
Only one response dummy predictor.
Chosen_none: grouped by chosen value
Lowhigh_none: grouped by 5 and 8, 12 and 25

Ver5_2:
6 binary predictors for risk and ambiguity level, do not separate for reward level, separate for each domain.
Only one response dummy predictor.

Ver5_3:
8 predictors for each domain: 4 outcome magnitude x 2 risk/ambig
16 predicotros altogether


Ver8:
Improve from Ver7. 
Only one response dummy predictor.
add_zero_columns_ver8


Ver8_1:
Only for Monetary domain
So that the empty columns for medical domain is removed
The existence of this version compromise the fact that Brainvoyager does not allow empty columns.
But Neuroelf is able to handle empty columns.


Ver9:
3 parametric modulators in the same model, with only mon and med for binary predicotr, plus one resposne binary predictor.
Parametric modulators include: 
Subjective rating
Ambiguity level
Risk level
Adjust_columns_ver9


