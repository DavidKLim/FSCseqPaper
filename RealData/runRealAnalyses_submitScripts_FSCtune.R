
# submit FSCtune runs as batch jobs into slurm

source("./RealData/runRealAnalyses_runScripts.R")
runScripts_FSC_tune(ncores=10, dataset="BRCA_pure", med_filt=500, MAD_filt=50, covariates='F', K_search=c(2:8), hrs=48)
runScripts_FSC_tune(ncores=10, dataset="BRCA_full", med_filt=500, MAD_filt=50, covariates='F', K_search=c(2:15), hrs=72)
runScripts_FSC_tune(ncores=15, dataset="BRCA_pure", med_filt=500, MAD_filt=50, covariates='T', K_search=c(2:8), hrs=72)
runScripts_FSC_tune(ncores=20, dataset="BRCA_full", med_filt=500, MAD_filt=50, covariates='T', K_search=c(2:15), hrs=96)
