The script in scripts/ runs GMM model fitting, best GMM selection, and filtering based on indiividual componenet AUC over a number of bins. 
Different AUC thresholds were tested. 
The script originally used summed AUCs for a model, almost no difference was found in the outcome for BOIs as seen in Table S1. 
The resulting peptidoform assignments by bin are split by AUC threshold in the out/ folder. 