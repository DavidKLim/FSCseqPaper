source("./Simulations/runSimAnalyses_createScripts.R")
source("./Simulations/runSimAnalyses_runScripts_slurm.R")

submit_jobs=function(method=c("FSC","iCl","MC","Others","NMF"),sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
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
            sim_index)
  }
  if("MC" %in% method){
    createScripts_MC(sigma_g,sigma_b, B, LFCb,
                      K, n, LFCg, pDEg, beta, phi,
                      K_search,sim_index)
    runMCs(sigma_g,sigma_b, B, LFCb,
           K, n, LFCg, pDEg, beta, phi,
           sim_index)
  }
  if("Others" %in% method){
    createScripts_Others(sigma_g,sigma_b, B, LFCb,
                         K, n, LFCg, pDEg, beta, phi,
                         K_search,sim_index)
    runOthers(sigma_g,sigma_b, B, LFCb,
               K, n, LFCg, pDEg, beta, phi,
               sim_index)
  }
  if("NMF" %in% method){
    createScripts_NMF(sigma_g,sigma_b, B, LFCb,
                         K, n, LFCg, pDEg, beta, phi,
                         K_search,sim_index)
    runNMFs(sigma_g,sigma_b, B, LFCb,
            K, n, LFCg, pDEg, beta, phi,
            sim_index)
  }
}

## changed existing directory names:
# mv nobatch B1_LFCb0
# mv batch1 B2_LFCb1
# mv batch2 B2_LFCb2
# mv batch3 B2_LFCb3

# submit all methods
submit_jobs(sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=2,n=c(50,100),sim_index=c(1:25))
submit_jobs(sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=4,n=c(100,200),sim_index=c(1:25))
# submit additional sim runs for main table
submit_jobs(method="FSC",sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=4,beta=12,sim_index=c(26:100))
# submit batch-affected simulation cases (LFCb=2,3), phi=0.35
submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=2,K=2,n=100,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=2,K=4,n=200,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=3,K=2,n=100,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=3,K=4,n=200,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
## phi=0.50
submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=2,K=2,n=100,LFCg=2,pDEg=0.05,beta=12,phi=0.50,sim_index=c(1:25))
submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=2,K=4,n=200,LFCg=2,pDEg=0.05,beta=12,phi=0.50,sim_index=c(1:25))
submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=3,K=2,n=100,LFCg=2,pDEg=0.05,beta=12,phi=0.50,sim_index=c(1:25))
submit_jobs(sigma_g=0.1,sigma_b=0,B=2,LFCb=3,K=4,n=200,LFCg=2,pDEg=0.05,beta=12,phi=0.50,sim_index=c(1:25))


########### EXTRA NMF RUNS (that were run separately) ############# (DONE)
# # submit all methods (DONE)    # do this first, then wait before doing batch
# submit_jobs(method="NMF",sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=2,n=c(50,100),sim_index=c(1:25))
# submit_jobs(method="NMF",sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=4,n=c(100,200),sim_index=c(1:25))
# # submit batch-affected simulation cases (LFCb=2,3) (DONE)
# submit_jobs(method="NMF",sigma_g=0.1,sigma_b=0,B=2,LFCb=2,K=2,n=100,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
# submit_jobs(method="NMF",sigma_g=0.1,sigma_b=0,B=2,LFCb=2,K=4,n=200,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
# submit_jobs(method="NMF",sigma_g=0.1,sigma_b=0,B=2,LFCb=3,K=2,n=100,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
# submit_jobs(method="NMF",sigma_g=0.1,sigma_b=0,B=2,LFCb=3,K=4,n=200,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))

########### phi=0.01 test (this is the case FSCseq did poorest in the main text, although better than other methods) ########### (SUBMITTED) #########
# simAllData function from SimulateRNASeqData.R script
simAllData(K=4,n_k=c(25,50),LFCg=c(1,2),pDEg=c(0.025,0.05),beta0=12,phi0=0.01,B=B,LFCb=LFCb,25)

submit_jobs(sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=4,n=100,LFCg=1,pDEg=0.025,beta=12,phi=0.01,sim_index=c(1:25))
submit_jobs(method=c("FSC","iCl","MC","Others","NMF"),sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=4,n=c(100,200),LFCg=c(1),pDEg=c(0.025,0.05),beta=12,phi=0.01,sim_index=c(1:25))

########### pDEg=0.50 test (this is the case FSCseq did poorest in the main text, although better than other methods) ########### (SUBMITTED) ##########
simAllData(K=4,n_k=c(25,50),LFCg=c(1,2),pDEg=0.5,beta0=12,phi0=c(0.15,0.35,0.5),B=B,LFCb=LFCb,25)

submit_jobs(sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=4,n=100,LFCg=1,pDEg=0.50,beta=12,phi=c(0.15,0.35,0.50),sim_index=c(1:25))
submit_jobs(method=c("FSC","iCl","MC","Others","NMF"),sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=4,n=c(100,200),LFCg=c(1),pDEg=0.5,beta=12,phi=c(0.15,0.35,0.50),sim_index=c(1:25))



# submit_jobs(method=c("iCl"),sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K=4,n=100,LFCg=1,pDEg=0.025,beta=12,phi=0.01,sim_index=c(1)) # seff 61142943

# # # tests NMF
# submit_jobs(method=c("NMF"),sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
#             K=c(2), n=c(50), LFCg=c(1), pDEg=c(0.025), beta=c(8), phi=c(0.15),
#             K_search='c(2:6)', lambda_search='seq(0.25,5,0.25)', alpha_search='c(0.01,seq(0.05,0.5,0.05))',
#             sim_index=c(2),ncores_FSC=1,ncores_iCl=5,nMB=5)      # seff 60659467: 40 mins, 980 MB         ####(INTNMF): 16.1 mins, 752 MB
# #
# submit_jobs(method=c("NMF"),sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
#             K=c(2), n=c(100), LFCg=c(1), pDEg=c(0.025), beta=c(8), phi=c(0.15),
#             K_search='c(2:6)', lambda_search='seq(0.25,5,0.25)', alpha_search='c(0.01,seq(0.05,0.5,0.05))',
#             sim_index=c(1),ncores_FSC=1,ncores_iCl=5,nMB=5)
# # # more 2_50_1.000000_0.025000_8.000000_0.150000_sim1_NMF.Rout   # seff 60659465: 126 mins, 1.05 GB                   #####(INTNMF):33.5 mins, 837 MB
# # # ll /pine/scr/d/e/deelim/out/Simulations/B1LFCb0/*_NMF.out
# #
# submit_jobs(method=c("NMF"),sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
#             K=c(4), n=c(200), LFCg=c(1), pDEg=c(0.025), beta=c(8), phi=c(0.15),
#             K_search='c(2:6)', lambda_search='seq(0.25,5,0.25)', alpha_search='c(0.01,seq(0.05,0.5,0.05))',
#             sim_index=c(1),ncores_FSC=1,ncores_iCl=5,nMB=5)     # seff 60659466: 151 mins, 1.05 GB                     ##### (INTNMF):39 mins, 891 MB
# #
#
# submit_jobs(method=c("NMF"),sigma_g=0.1,sigma_b=0, B=2, LFCb=2,
#             K=c(2), n=c(100), LFCg=c(2), pDEg=c(0.05), beta=c(12), phi=c(0.35),
#             K_search='c(2:6)',sim_index=c(1))      # seff 60659484 38 mins, 1.01 GB

## running interactively
# library(IntNMF)
# source("runSimAnalyses.R")
# runSimAnalyses_NMF(K=2,n=50,LFCg=1,pDEg=0.025,beta=8,phi=0.15,
#                    sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K_search=c(2:6),sim_index=1)
# system("ls /pine/scr/d/e/deelim/out/Simulations/B1_LFCb0/*_NMF.out")
# load("/pine/scr/d/e/deelim/out/Simulations/B1_LFCb0/2_50_1.000000_0.025000_8.000000_0.150000_sim1_NMF.out")
# ls()
# NMF_summary$cls_metrics

# ## running interactively
# source("runSimAnalyses.R")
# n=c(100,200);LFCg=c(1,2);pDEg=c(0.025,0.05);phi=0.01;sim_index=c(1:25)
# for(a in 1:length(n)){for(b in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(phi)){for(f in 1:length(sim_index)){
#   if(file.exists(sprintf("/pine/scr/d/e/deelim/out/Simulations/0.100000_0.000000/B1_LFCb0/4_%d_%f_%f_12.000000_%f_sim%d_res_FSCpred.out",
#                          n[a],LFCg[b],pDEg[d],phi[e],sim_index[f]))){ next }
#   runSimAnalyses_FSC(K=4,n=n[a],LFCg=LFCg[b],pDEg=pDEg[d],beta=12,phi=phi[e],
#                      sigma_g=0.1,sigma_b=0,B=1,LFCb=0,
#                      K_search=c(2:6),lambda_search=seq(0.25,5,0.25),alpha_search=c(0.01,seq(0.05,0.5,0.05)),
#                      sim_index=sim_index[f],ncores=1,nMB=5,trace=T)
# }}}}}
# n=c(100,200);LFCg=c(1,2);pDEg=c(0.5);phi=c(0.15,0.35,0.5);sim_index=c(1:25)
# for(a in 1:length(n)){for(b in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(phi)){for(f in 1:length(sim_index)){
#   if(file.exists(sprintf("/pine/scr/d/e/deelim/out/Simulations/0.100000_0.000000/B1_LFCb0/4_%d_%f_%f_12.000000_%f_sim%d_res_FSCpred.out",
#                          n[a],LFCg[b],pDEg[d],phi[e],sim_index[f]))){ next }
#   runSimAnalyses_FSC(K=4,n=n[a],LFCg=LFCg[b],pDEg=pDEg[d],beta=12,phi=phi[e],
#                      sigma_g=0.1,sigma_b=0,B=1,LFCb=0,
#                      K_search=c(2:6),lambda_search=seq(0.25,5,0.25),alpha_search=c(0.01,seq(0.05,0.5,0.05)),
#                      sim_index=sim_index[f],ncores=1,nMB=5,trace=T)
# }}}}}
#
# runSimAnalyses_FSC(K=4,n=200,LFCg=1,pDEg=0.05,beta=12,phi=0.01,
#                    sigma_g=0.1,sigma_b=0,B=1,LFCb=0,
#                    K_search=c(2:6),lambda_search=seq(0.25,5,0.25),alpha_search=c(0.01,seq(0.05,0.5,0.05)),
#                    sim_index=1,ncores=1,nMB=5,trace=T)
# load("/pine/scr/d/e/deelim/out/Simulations/0.100000_0.000000/B1_LFCb0/4_200_1.000000_0.050000_12.000000_0.010000_sim1_res_FSCpred.out")
# ls()
# NMF_summary$cls_metrics

