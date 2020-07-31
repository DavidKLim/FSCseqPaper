# source functions that perform each comparative method
# source("BatchPred/runSimAnalyses_functions.R")

## Wrapper for FSCseq
runSimAnalyses_FSC = function(K=2,n=100,LFCg=1,pDEg=0.05,beta=8,phi=0.15,
                              sigma_g=0.1,sigma_b=0,B=1,LFCb=0,
                              K_search=c(2:6),lambda_search=seq(0.25,5,0.25),alpha_search=c(0.01,seq(0.05,0.5,0.05)),
                              sim_index=1,ncores=1,nMB=5,trace=T){
  library(fpc)
  library(mclust)
  library(cluster)
  library(NbClust)
  library(FSCseq)
  library(iClusterPlus)
  library(NB.MClust)
  library(pryr)
  library(SummarizedExperiment)

  dir_name=sprintf("./Simulations/%f_%f/B%d_LFCb%d/BatchPred",sigma_g,sigma_b,B,LFCb)

  # use file.name to specify data.file to load, and trace.folder to create to store diagnostics
  sim.cond = sprintf("%d_%d_%f_%f_%f_%f_sim%d",K,n,LFCg,pDEg,beta,phi,sim_index)

  file.name = sprintf("%s/%s",dir_name,sim.cond)
  data.file = sprintf("%s_data.RData",file.name)

  trace.folder = sprintf("%s/%s",dir_name,sim.cond)
  ifelse(!dir.exists(trace.folder),
         dir.create(trace.folder),
         FALSE)

  # if sim has already been run (results exist in Simulations folder) --> don't run
  res.file = sprintf("%s/%s_res",dir_name,sim.cond)

  load(data.file)
  # Input simulated data
  cts = sim.dat$cts
  cls = sim.dat$cls
  batch = sim.dat$batch
  B= sim.dat$sim_params$B
  g= sim.dat$sim_params$g
  DEg_ID = sim.dat$DEg_ID

  res1=sprintf("%s_FSC.out",res.file)
  if(!file.exists(res1)){
    cat("Running FSCseq analysis...\n")
    FSC_summary = run_FSCseq(sim.dat=sim.dat, batch=batch, K_search=K_search,
                             lambda_search=lambda_search, alpha_search=alpha_search,
                             trace.folder=trace.folder, trace=trace,ncores=ncores, nMB=nMB)
    save(FSC_summary,file=res1)
  } else {
    cat("\nFSCseq analysis already run...\n")
    load(res1)
  }

  #### PREDICTION ON SIMULATED DATASET ####
  res7=sprintf("%s_FSCpred.out",res.file)
  if(!file.exists(res7)){
    cat("\nAnalyzing simulated prediction data via FSCseq...\n")
    FSC_predict_summary = run_FSCseq_predict(sim.dat=sim.dat,FSC_summary=FSC_summary)
    save(FSC_predict_summary,file=res7)
  } else{
    cat("\nFSCpred already run...\n")
  }

  # plot heatmaps
  norm_y=FSC_summary$norm_y; idx=FSC_summary$idx
  library(pheatmap)
  norm_y2=norm_y[idx,]; DEg_ID2=DEg_ID[idx]
  colnames(norm_y2)=c(1:ncol(norm_y2))
  rownames(norm_y2)=c(1:nrow(norm_y2))
  cls2=as.numeric(as.factor(FSC_summary$cls))
  annotation_col = data.frame(factor(cls),
                              factor(cls2))
  colnames(annotation_col)=c("truth","FSC")
  rownames(annotation_col)=colnames(norm_y2)
  newCols <- colorRampPalette(grDevices::rainbow(max(length(unique(cls)),length(unique(cls2)))))
  mycolors_true=newCols(length(unique(cls))); names(mycolors_true)=unique(cls)[order(unique(cls))]
  mycolors_FSC=newCols(length(unique(cls2))); names(mycolors_FSC)=unique(cls2)[order(unique(cls2))]
  mycolors=list(truth=mycolors_true,FSC=mycolors_FSC)

  annotation_row=data.frame(as.factor(abs(FSC_summary$fit$phi-phi)>0.1))
  colnames(annotation_row)="Overestimated phi"
  rownames(annotation_row)=rownames(norm_y2)
  mycolors_1=c("#FFFFFF","#000000"); names(mycolors_1)=c("TRUE","FALSE")
  mycolors2=list("Overestimated phi"=mycolors_1)
  png(sprintf("%s/HM_disc.png",trace.folder),height=600,width=600)
  pheatmap(log(norm_y2[DEg_ID2,order(cls,cls2)]+0.1),cluster_cols=F,annotation_col=annotation_col,annotation_row=annotation_row,annotation_colors=c(mycolors,mycolors2),main="FSC disc genes")
  dev.off()
  png(sprintf("%s/HM_nondisc.png",trace.folder),height=600,width=600)
  pheatmap(log(norm_y2[!DEg_ID2,order(cls,cls2)][1:500,]+0.1),cluster_cols=F,annotation_col=annotation_col,annotation_row=annotation_row,annotation_colors=c(mycolors,mycolors2),main="FSC 500 random nondisc genes")
  dev.off()
  png(sprintf("%s/HM_coefs_disc.png",trace.folder),height=600,width=600)
  pheatmap(FSC_summary$fit$coefs[DEg_ID2,],cluster_cols=F,cluster_rows=F,main="FSC coefficients, disc genes")
  dev.off()
  png(sprintf("%s/HM_coefs_nondisc.png",trace.folder),height=600,width=600)
  pheatmap(FSC_summary$fit$coefs[!DEg_ID2,],cluster_cols=F,cluster_rows=F,main="FSC coefficients, nondisc genes")
  dev.off()
  png(sprintf("%s/phi_ests.png",trace.folder),height=600,width=600)
  plot(FSC_summary$fit$phi,main=sprintf("Phi ests vs. gene index. True phi=%f",phi))
  dev.off()

}

# Implementation of FSCseq (READY)
run_FSCseq=function(sim.dat,batch,K_search,
                    lambda_search,alpha_search,
                    trace.folder,ncores,nMB,trace){

  cts=sim.dat$cts; cls=sim.dat$cls; DEg_ID=sim.dat$DEg_ID
  library(FSCseq)
  start_FSC = Sys.time()
  n_rinits=1; method="CEM"
  FSCseq_results = FSCseq::FSCseq_workflow(cts=cts,ncores=ncores,batch=batch,true_cls=cls,true_disc=DEg_ID,
                                           method="CEM",n_rinits=1,med_filt=100,MAD_filt=50,
                                           K_search=K_search,lambda_search=lambda_search,alpha_search=alpha_search,
                                           OS_save=TRUE,trace=trace,trace.prefix="",nMB=nMB,dir_name=trace.folder)

  cat("Tuning complete.\n")
  # save BICs, n_its, and time_elap from tuning
  # tuning_results = list(BICs=BICs,num_est_params=num_est_params,
  #                       n_its=n_its,
  #                       time_elap=time_elap,
  #                       cls=tuning_cls,
  #                       TPRs=tuning_TPRs,
  #                       FPRs=tuning_FPRs,
  #                       Ks=tuning_K,
  #                       ASWs=ASWs)

  fit=FSCseq_results$results$fit
  norm_y=FSCseq_results$processed.dat$norm_y; idx=FSCseq_results$processed.dat$idx

  max_lambda=fit$lambda; max_alpha=fit$alpha
  end_FSC = Sys.time()

  cat(paste("Final fit complete. Total time elapsed: ",as.numeric(end_FSC-start_FSC,units="secs"),"s.\n",sep=""))

  # Calculate sensitivity (TPR) and false positive rate (FPR)
  FSC_discriminatory=fit$discriminatory
  TPR = sum(DEg_ID[idx] & FSC_discriminatory)/sum(DEg_ID[idx])     # true disc found to truly be disc
  FPR = sum(!DEg_ID[idx] & FSC_discriminatory)/sum(!DEg_ID[idx])   # not disc found falsely to be disc

  # Calculate cls metrics
  cls_metrics = cluster.stats(d=NULL,fit$clusters,cls,compareonly=T)

  return(list(fit=fit, processed.dat=processed.dat, results = FSCseq_results$results,   ## processed.dat and results aren't required, but works better rn with FSCseq_predict_workflow
              K=length(unique(fit$clusters)),
              max_lambda=max_lambda,max_alpha=max_alpha,
              cls=fit$clusters,
              cls_metrics=cls_metrics,
              time_elap=as.numeric(end_FSC-start_FSC,units="secs"),
              TPR=TPR,FPR=FPR,
              norm_y=norm_y,idx=idx))
}
# Implementation of FSCseq-prediction (READY)
run_FSCseq_predict=function(sim.dat, FSC_summary){
  library(FSCseq)
  cts_train=sim.dat$cts; cts_pred=sim.dat$cts_pred; cls_pred=sim.dat$cls_pred; batch_train=sim.dat$batch; batch_pred=sim.dat$batch_pred

  FSCseq_predict_results = FSCseq::FSCseq_predict_workflow(res=FSC_summary, X_covar_train = NULL, cts_train=cts_train, SF_train=NULL, batch_train=batch_train,
                                                           X_covar_pred = NULL, cts_pred=cts_pred, batch_pred=batch_pred)

  res = FSCseq_predict_results
  pARI = adjustedRandIndex(res$pred_cls,cls_pred)

  cls_metrics=cluster.stats(d=NULL,res$pred_cls,cls_pred,compareonly=T)
  return(list(fit=res,
              cls=res$clusters,
              pARI=pARI,
              cls_metrics=cls_metrics))
}

## Function to create the relevant scripts
createScripts_FSC=function(sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                 K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                 K_search='c(2:6)', lambda_search='seq(0.25,5,0.25)', alpha_search='c(0.01,seq(0.05,0.5,0.05))',
                 sim_index=c(1:25),ncores=1,nMB=5){

  if(B==1){LFCb=0}

  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d/BatchPred",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d/BatchPred",sigma_g,sigma_b,B,LFCb)
  dir.create(dir_name1, recursive=TRUE, showWarnings = F)
  dir.create(dir_name2, recursive=TRUE, showWarnings = F)

  for(a in 1:length(K)){for(b in 1:length(n)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta)){for(f in 1:length(phi)){for(g in 1:length(sim_index)){
    fname = sprintf("%d_%d_%f_%f_%f_%f_sim%d",
                    K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],
                    sim_index[g])

    out = sprintf("%s/%s",dir_name1,fname)
    out_script = sprintf("%s/%s_FSC",dir_name2,fname)

    res1=sprintf("%s_res_FSC.out",out)
    res2=sprintf("%s_res_FSCpred.out",out)

    # submit job only if the FSC .out files do not exist
    if(!file.exists(res1) & !file.exists(res2)){
      cmd = rep(0, 2)
      cmd[1] = "unlink('.RData') \n source('runSimAnalyses_BatchPred.R') \n"
      cmd[2] = sprintf("res = runSimAnalyses_FSC(K=%d,n=%d,LFCg=%f,pDEg=%f,beta=%f,phi=%f,
                            sigma_g=%f,sigma_b=%f,B=%d,LFCb=%d,
                            K_search=%s,lambda_search=%s,alpha_search=%s,
                            sim_index=%d,ncores=%d,nMB=%d)\n",
                       K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],
                       sigma_g,sigma_b,B,LFCb,
                       K_search,lambda_search,alpha_search,
                       sim_index[g],ncores,nMB)
      cmdf = paste(cmd, collapse = "")
      write.table(cmdf, file = out_script, col.names = F, row.names = F, quote = F)
    }
  }}}}}}}

}

## Function to submit jobs to slurm
runFSCs=function(ncores=1,sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                 K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                 sim_index=c(1:25)){

  if(B==1){LFCb=0}
  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d/BatchPred",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d/BatchPred",sigma_g,sigma_b,B,LFCb)

  for(a in 1:length(K)){for(b in 1:length(n)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta)){for(f in 1:length(phi)){for(g in 1:length(sim_index)){
    fname = sprintf("%d_%d_%f_%f_%f_%f_sim%d",
                    K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],
                    sim_index[g])

    out = sprintf("%s/%s",dir_name1,fname)
    out_script = sprintf("%s/%s_FSC",dir_name2,fname)

    res1=sprintf("%s_res_FSC.out",out)
    res2=sprintf("%s_res_FSCpred.out",out)

    # submit job only if the FSC .out files do not exist
    if(!file.exists(res1) & !file.exists(res2)){

      mem = 10000 + 500*(n[b]==200)^2
      hrs = 24*(n[b]==50) + 36*(n[b]==100) + 48*(n[b]==200)

      run = sprintf("sbatch -p general -N 1 -n %d --mem=%d -t %d:00:00 -J %s_FSC --wrap='R CMD BATCH %s'",ncores,mem,hrs,fname,out_script)
      print(run)
      Sys.sleep(0.2)
      system(run)
    }

  }}}}}}}
}

## Function to create scripts (createScripts), then submit the jobs (runFSCs) (READY)
submit_jobs=function(method=c("FSC"),sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                     K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                     K_search='c(2:6)', lambda_search='seq(0.25,5,0.25)', alpha_search='c(0.01,seq(0.05,0.5,0.05))',
                     sim_index=c(1:25),ncores_FSC=1,nMB=5 #,ncores_iCl=5
){
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
  # if("iCl" %in% method){
  #   createScripts_iCl(ncores_iCl,sigma_g,sigma_b, B, LFCb,
  #                     K, n, LFCg, pDEg, beta, phi,
  #                     K_search,sim_index)
  #   runiCls(ncores_iCl,sigma_g,sigma_b, B, LFCb,
  #           K, n, LFCg, pDEg, beta, phi,
  #           sim_index)
  # }
  # if("MC" %in% method){
  #   createScripts_MC(sigma_g,sigma_b, B, LFCb,
  #                    K, n, LFCg, pDEg, beta, phi,
  #                    K_search,sim_index)
  #   runMCs(sigma_g,sigma_b, B, LFCb,
  #          K, n, LFCg, pDEg, beta, phi,
  #          sim_index)
  # }
  # if("Others" %in% method){
  #   createScripts_Others(sigma_g,sigma_b, B, LFCb,
  #                        K, n, LFCg, pDEg, beta, phi,
  #                        K_search,sim_index)
  #   runOthers(sigma_g,sigma_b, B, LFCb,
  #             K, n, LFCg, pDEg, beta, phi,
  #             sim_index)
  # }
  # if("NMF" %in% method){
  #   createScripts_NMF(sigma_g,sigma_b, B, LFCb,
  #                     K, n, LFCg, pDEg, beta, phi,
  #                     K_search,sim_index)
  #   runNMFs(sigma_g,sigma_b, B, LFCb,
  #           K, n, LFCg, pDEg, beta, phi,
  #           sim_index)
  # }
}
###########################################################################
# Simulate the data
## run SimulateData/BatchPred/SimulateRNASeqData_BatchPred.R

# Submit the jobs (SUBMITTED)
source("runSimAnalyses_BatchPred.R")
submit_jobs(method=c("FSC"),sigma_g=0.1,sigma_b=0,B=2,LFCb=2,K=2,n=100,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
submit_jobs(method=c("FSC"),sigma_g=0.1,sigma_b=0,B=2,LFCb=2,K=4,n=200,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
submit_jobs(method=c("FSC"),sigma_g=0.1,sigma_b=0,B=2,LFCb=3,K=2,n=100,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))
submit_jobs(method=c("FSC"),sigma_g=0.1,sigma_b=0,B=2,LFCb=3,K=4,n=200,LFCg=2,pDEg=0.05,beta=12,phi=0.35,sim_index=c(1:25))

source("runSimAnalyses_BatchPred.R")
library(mclust); library(fpc)
for(s in c(1:25)){
  for(LFCb in c(2,3)){for(K in c(2,4)){
  B=2; n=50*K
  res.file = sprintf("Simulations/0.100000_0.000000/B%d_LFCb%d/BatchPred/%d_%d_2.000000_0.050000_12.000000_0.350000_sim%d",
                     B,LFCb,K,n,s)
  if(file.exists(sprintf("%s_res_FSC.out",res.file))){
    load(sprintf("%s_res_FSC.out",res.file))
    load(sprintf("%s_data.RData",res.file))
    print(res.file)
    cat("\nAnalyzing simulated prediction data via FSCseq...\n")
    FSC_predict_summary = run_FSCseq_predict(sim.dat=sim.dat,FSC_summary=FSC_summary)
    save(FSC_predict_summary,file=sprintf("%s_res_FSCpred.out",res.file))
    print(paste("ARI (train): ", adjustedRandIndex(FSC_summary$cls, sim.dat$cls)))
    print(paste("ARI (pred): ", adjustedRandIndex(FSC_predict_summary$fit$pred_cls, sim.dat$cls_pred)))
  }
  }}
}




########## COLLECTSIMS ########
collectSims_BatchPred=function(sigma_g=0.1,sigma_b=0,B=1,LFCb=0,
                            K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                            compared_methods=c("FSC","pFSC","iCl","HC","KM","NBMB","NMF","lMC","vMC","rMC"),
                            sim_index=c(1:25)){

  dir_name=sprintf("./Simulations/%f_%f/B%d_LFCb%d/BatchPred",sigma_g,sigma_b,B,LFCb)
  library(mclust)

  metrics=c("K","OA","time","ARI")
  metrics_pred=c("ARI")

  n_FSC_disc = 2 * ("FSC" %in% compared_methods)
  n_metrics = length(compared_methods[compared_methods!="pFSC"])*length(metrics)
  n_metrics_pred = length(metrics_pred) * ("pFSC" %in% compared_methods)
  n_penalty = 2 * ("FSC" %in% compared_methods)
  n_iCl_disc = 2 * ("iCl" %in% compared_methods)
  # data.frame to store all results
  n_conditions = length(K)*length(n)*length(LFCg)*length(pDEg)*length(beta)*length(phi)
  df = data.frame(matrix(ncol = (6 + n_FSC_disc + n_metrics + n_metrics_pred + n_penalty + n_iCl_disc), nrow = n_conditions))

  # 6 variable sim params (truth)
  colnames(df)[1:6]=c("K","n","LFCg","pDEg","beta","phi")

  # 2: TPR, FPR for FSCseq feature selection
  index=6+1
  if(n_FSC_disc==2){colnames(df)[c(index,index+1)] = c("TPR_FSC","FPR_FSC")}

  # 32: (K, OA, time, and ARI) for all 8 methods
  index=6 + n_FSC_disc + 1
  for(j in 1:length(compared_methods)){
    if(compared_methods[j]=="pFSC"){
      for(i in 1:length(metrics_pred)){
        colnames(df)[index] = sprintf("%s_%s",metrics_pred[i],compared_methods[j])
        index=index+1
      }
    }else{
      for(i in 1:length(metrics)){
        colnames(df)[index] = sprintf("%s_%s",metrics[i],compared_methods[j])
        index=index+1
      }
    }
  }

  index=6 + n_FSC_disc + n_metrics + n_metrics_pred + 1
  if(n_penalty==2){colnames(df)[c(index,index+1)] = c("lambda","alpha")}
  index=6 + n_FSC_disc + n_metrics + n_metrics_pred + n_penalty + 1
  if(n_iCl_disc==2){colnames(df)[c(index,index+1)] = c("TPR_iCl","FPR_iCl")}

  summary_methods = paste(compared_methods,"_summary",sep="")
  calculate_metrics = function(metric,summary_obj,cls,true.K,batch=NULL){
    summary=eval(parse(text=summary_obj))
    if(metric=="K"){
      result=length(unique(summary$cls))
    } else if(metric=="OA"){
      result=((length(unique(summary$cls))==true.K)^2)
    } else if(metric=="time"){
      result=(summary$time_elap)
    } else if(metric=="ARI"){
      result=(summary$cls_metrics$corrected.rand)
    }


    if(is.null(result)){return(NA)} else{return(result)}
  }

  # index for each combination of sim params
  index=1
  all_sims_df = data.frame(matrix(nrow=0,ncol=ncol(df)))


  for(a in 1:length(K)){for(b in 1:length(n)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta)){for(f in 1:length(phi)){
    # tabulate all results for all sim datasets for each condition
    all_sims_mat = data.frame(matrix(NA,nrow=length(sim_index),ncol=(ncol(df))))
    colnames(all_sims_mat)=colnames(df)
    for(g in 1:length(sim_index)){
      # read in main res out file and MC_summary out file
      fname = sprintf("%d_%d_%f_%f_%f_%f_sim%d_res",
                      K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],
                      sim_index[g])
      res.file=sprintf("%s/%d_%d_%f_%f_%f_%f_sim%d",
                       dir_name,K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],sim_index[g])
      res1=sprintf("%s_res_FSC.out",res.file)
      if(file.exists(res1) & ("FSC" %in% compared_methods)){load(res1)}else if(!file.exists(res1) & ("FSC" %in% compared_methods)){next}
      res7=sprintf("%s_res_FSCpred.out",res.file)
      if(file.exists(res7) & ("pFSC" %in% compared_methods)){
        load(res7)
        pFSC_summary=FSC_predict_summary
      }else if(!file.exists(res7) & ("pFSC" %in% compared_methods)){next}
      res4=sprintf("%s_KM.out",res.file)
      if(file.exists(res4) & ("KM" %in% compared_methods)){load(res4)}else if(!file.exists(res4) & ("KM" %in% compared_methods)){next}
      res5=sprintf("%s_NBMB.out",res.file)
      if(file.exists(res5) & ("NBMB" %in% compared_methods)){load(res5)}else if(!file.exists(res5) & ("NBMB" %in% compared_methods)){next}
      res3=sprintf("%s_HC.out",res.file)
      if(file.exists(res3) & ("HC" %in% compared_methods)){load(res3)}else if(!file.exists(res3) & ("HC" %in% compared_methods)){next}
      res2=sprintf("%s_iCl.out",res.file)
      if(file.exists(res2) & ("iCl" %in% compared_methods)){load(res2)}else if(!file.exists(res2) & ("iCl" %in% compared_methods)){next}
      res6=sprintf("%s_MC.out",res.file)
      if(file.exists(res6) & (all(c("lMC","vMC","rMC") %in% compared_methods))){
        load(res6)
        lMC_summary=MC_summary$lMC; vMC_summary=MC_summary$vMC; rMC_summary=MC_summary$rMC
      }else if(!file.exists(res6) & (all(c("lMC","vMC","rMC") %in% compared_methods))){next}

      res8=sprintf("%s_NMF.out",res.file)
      if(file.exists(res8) & ("NMF" %in% compared_methods)){load(res8)}else if(!file.exists(res8) & ("NMF" %in% compared_methods)){next}

      print(res.file)

      load(sprintf("%s/%d_%d_%f_%f_%f_%f_sim%d_data.RData",
                   dir_name,K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],sim_index[g]))
      cls = sim.dat$cls
      batch = sim.dat$batch
      cls_pred=sim.dat$cls_pred

      # iterate over all sim runs with same condition, create huge matrix

      all_sims_mat[g,1:6] = c(K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f])

      index2=6+1
      # TPR and FPR (just FSC)
      if("FSC" %in% compared_methods){all_sims_mat[g,c(index2,index2+1)] = c(FSC_summary$TPR,FSC_summary$FPR)}

      # index2 for each sim dataset within fixed combination of sim params
      index2=6 + n_FSC_disc + 1

      for(j in 1:length(summary_methods)){
        if(summary_methods[j]=="pFSC_summary"){
          for(i in 1:length(metrics_pred)){
            all_sims_mat[g,index2] = calculate_metrics(metrics_pred[i],summary_methods[j],cls_pred,K[a],batch)
            index2=index2+1
          }
        } else{
          for(i in 1:length(metrics)){
            all_sims_mat[g,index2] = calculate_metrics(metrics[i],summary_methods[j],cls,K[a],batch)
            index2=index2+1
          }
        }
      }
      index2=6 + n_FSC_disc + n_metrics + n_metrics_pred + 1
      if("FSC" %in% compared_methods){all_sims_mat[g,c(index2,index2+1)] = c(FSC_summary$fit$lambda,FSC_summary$fit$alpha)}
      index2=6 + n_FSC_disc + n_metrics + n_metrics_pred + n_penalty + 1
      if("iCl" %in% compared_methods){all_sims_mat[g,c(index2,index2+1)] = c(iCl_summary$TPR,iCl_summary$FPR)}

    }
    df[index,1:ncol(df)]=colMeans(all_sims_mat, na.rm=T)
    all_sims_df = rbind(all_sims_df,all_sims_mat)

    # average out every variable
    index = index+1
  }}}}}}


  return(list(df=df,all_sims_df=all_sims_df
  )) # return result to evaluate interactively
}

res.dir="./Simulations/Results"
ifelse(!dir.exists(res.dir),dir.create(res.dir),FALSE)

df2batch2=collectSims_BatchPred(sigma_g=0.1,sigma_b=0,B=2,LFCb=2,
                      K=c(2), n=c(100), LFCg=c(2), pDEg=c(0.05), beta=c(12), phi=c(0.35),
                      compared_methods=c("FSC",'pFSC'),
                      sim_index=c(1:25))
df2batch2$df$gamma=2;df2batch2$all_sims_df$gamma=2
df4batch2=collectSims_BatchPred(sigma_g=0.1,sigma_b=0,B=2,LFCb=2,
                      K=c(4), n=c(200), LFCg=c(2), pDEg=c(0.05), beta=c(12), phi=c(0.35),
                      compared_methods=c("FSC",'pFSC'),
                      sim_index=c(1:25))
df4batch2$df$gamma=2;df4batch2$all_sims_df$gamma=2

df2batch3=collectSims_BatchPred(sigma_g=0.1,sigma_b=0,B=2,LFCb=3,
                      K=c(2), n=c(100), LFCg=c(2), pDEg=c(0.05), beta=c(12), phi=c(0.35),
                      compared_methods=c("FSC",'pFSC'),
                      sim_index=c(1:25))
df2batch3$df$gamma=3;df2batch3$all_sims_df$gamma=3
df4batch3=collectSims_BatchPred(sigma_g=0.1,sigma_b=0,B=2,LFCb=3,
                      K=c(4), n=c(200), LFCg=c(2), pDEg=c(0.05), beta=c(12), phi=c(0.35),
                      compared_methods=c("FSC",'pFSC'),
                      sim_index=c(1:25))
df4batch3$df$gamma=3;df4batch3$all_sims_df$gamma=3
df2batch=list(df=rbind(df2batch2$df,df2batch3$df),all_sims_df=rbind(df2batch2$all_sims_df,df2batch3$all_sims_df))
df4batch=list(df=rbind(df4batch2$df,df4batch3$df),all_sims_df=rbind(df4batch2$all_sims_df,df4batch3$all_sims_df))

save(df2batch,file=sprintf("%s/df2batch_BatchPred.out",res.dir)); save(df4batch,file=sprintf("%s/df4batch_BatchPred.out",res.dir))
