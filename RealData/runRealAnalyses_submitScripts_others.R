
# submit FSCtune runs as batch jobs into slurm

source("./RealData/runRealAnalyses_runScripts.R")
runScripts_iCl(ncores=13, dataset="BRCA_pure", med_filt=500, MAD_filt=50, K_search=c(2:8), hrs=48)
runScripts_iCl(ncores=25, dataset="BRCA_full", med_filt=500, MAD_filt=50, K_search=c(2:15), hrs=48)

runScripts_MC(dataset="BRCA_pure", med_filt=500, MAD_filt=50, K_search=c(2:8), hrs=24)
runScripts_MC(dataset="BRCA_full", med_filt=500, MAD_filt=50, K_search=c(2:15), hrs=48)

runScripts_Others(dataset="BRCA_pure", med_filt=500, MAD_filt=50, K_search=c(2:8), hrs=24)
runScripts_Others(dataset="BRCA_full", med_filt=500, MAD_filt=50, K_search=c(2:15), hrs=48)
