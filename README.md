# BayesFactorsPsych
 
Materials for a comment submitted to ` `.

The main script is case1-fulldata.R It uses the case 1 data (full vs aggregate data) to explore the use of LOO-CV in various forms with data simulated using the target article's method (in target-article-code).

The results are saved in the results folder, which contains:

- kfold-grp.rds (kfold leaving each respondent out)
- kfold-grp-correctse (kfold with se's adjusted as per https://avehtari.github.io/modelselection/rats_kcv.html#53_Grouped_K-fold_for_leave-one-group-out
)
- kfold-ind.rds (traditional leave one obs out, approximated using LOO-PSIS)
- model-fits-fulldata.rds (model fits for models 3 - 6 using the full data)