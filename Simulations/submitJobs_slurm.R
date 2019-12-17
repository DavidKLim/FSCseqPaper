source("./Simulations/runSimAnalyses_createScripts.R")
source("./Simulations/runSimAnalyses_runScripts_slurm.R")

submit_jobs=function(method=c("FSC","iCl","MC","Others"),sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
            K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
            K_search='c(2:6)', lambda_search='seq(0.25,5,0.25)', alpha_search='c(0.01,seq(0.05,0.5,0.05))',
            sim_index=c(1:25),ncores_FSC=1,ncores_iCl=5,nMB=5){
  # method="Others" includes HC, KM, and NBMB

  if("FSC" %in% method){
    createScripts_FSC(sigma_g,sigma_b, B, LFCb,
                      K, n, LFCg, pDEg, beta, phi,
                      K_search, lambda_search, alpha_search,
                      sim_index,ncores_FSC,nMB)
    runFSCs(ncores_FSC,sigma_g,sigma_b, B, LFCb,
            K, n, LFCg, pDEg, beta, phi,
            sim_index)
  }
  if("iCl" %in% method){
    createScripts_iCl(ncores_iCl,sigma_g,sigma_b, B, LFCb,
                       K, n, LFCg, pDEg, beta, phi,
                       K_search,sim_index)
    runiCls(ncores_iCl,sigma_g,sigma_b, B, LFCb,
            K, n, LFCg, pDEg, beta, phi,
            K_search,sim_index)
  }
  if("MC" %in% method){
    createScripts_MC(sigma_g,sigma_b, B, LFCb,
                      K, n, LFCg, pDEg, beta, phi,
                      K_search,sim_index)
    runMCs(sigma_g,sigma_b, B, LFCb,
           K, n, LFCg, pDEg, beta, phi,
           K_search,sim_index)
  }
  if("Others" %in% method){
    createScripts_Others(sigma_g,sigma_b, B, LFCb,
                         K, n, LFCg, pDEg, beta, phi,
                         K_search,sim_index)
    runOthers(sigma_g,sigma_b, B, LFCb,
               K, n, LFCg, pDEg, beta, phi,
               K_search,sim_index)
  }
}

# submit all methods
submit_jobs(sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=2,n=c(50,100),sim_index=c(1:25))
submit_jobs(sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=4,n=c(100,200),sim_index=c(1:25))

# submit additional sim runs for main table
submit_jobs(method="FSC",sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=4,beta=12,sim_index=c(26:100))

# submit batch-affected simulation cases (LFCb=2,3)
submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=2,K=2,n=100,LFCg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=2,K=4,n=200,LFCg=0.05,beta=12,phi=0.35,sim_index=c(1:25))

submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=3,K=2,n=100,LFCg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=3,K=4,n=200,LFCg=0.05,beta=12,phi=0.35,sim_index=c(1:25))

