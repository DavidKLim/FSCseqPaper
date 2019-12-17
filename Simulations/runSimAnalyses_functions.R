# Functions that perform each comparative analysis. Sourced in wrapper runSimAnalyses.R

## NEED TO EDIT FSCseq TO INCORPORATE PACKAGE
run_FSCseq=function(sim.dat,batch,K_search,
                    lambda_search,alpha_search,
                    trace.folder,ncores,nMB,trace){

  cts=sim.dat$cts; cls=sim.dat$cls; DEg_ID=sim.dat$DEg_ID
  library(FSCseq)
  start_FSC = as.numeric(Sys.time())
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
  end_FSC = as.numeric(Sys.time())

  cat(paste("Final fit complete. Total time elapsed: ",end_FSC-start_FSC,"s.\n",sep=""))

  # Calculate sensitivity (TPR) and false positive rate (FPR)
  FSC_discriminatory=fit$discriminatory
  TPR = sum(DEg_ID & FSC_discriminatory)/sum(DEg_ID)     # true disc found to truly be disc
  FPR = sum(!DEg_ID & FSC_discriminatory)/sum(!DEg_ID)   # not disc found falsely to be disc

  # Calculate cls metrics
  cls_metrics = cluster.stats(d=NULL,fit$clusters,cls,compareonly=T)

  return(list(fit=fit,
              K=length(unique(fit$clusters)),
              max_lambda=max_lambda,max_alpha=max_alpha,
              cls=fit$clusters,
              cls_metrics=cls_metrics,
              time_elap=(end_FSC-start_FSC),
              TPR=TPR,FPR=FPR,
              norm_y=norm_y,idx=idx))
}
run_iCl=function(cts,true_cls,K_search,
                 ncores,n.lambda,DEg_ID){

  start_iCl = as.numeric(Sys.time())
  if(K_search[1]==1){K_search=K_search[-1]}  # iCl can't compare K=1 --> remove this

  iCl_K_search = K_search-1      # iCl_K = K - 1 --> K = iCl_K + 1

  cts=round(cts,0)       # coerce normalized counts to integer

  # Tuning
  iClust_OS <- list()
  for(c in 1:length(iCl_K_search)){
    cv.fit = tune.iClusterPlus(cpus=ncores,dt1=t(cts),type="poisson",K=iCl_K_search[c],n.lambda=n.lambda,maxiter=20)
    iClust_OS[[c]] = cv.fit
  }
  BIC_mat = getBIC(iClust_OS)
  dev_mat = getDevR(iClust_OS)
  lambda_ids = apply(BIC_mat,2,which.min)   # indices of optimal lambda for each number of clusters
  devs=rep(0,length(iCl_K_search))
  for(c in 1:length(iCl_K_search)){
    devs[c] = dev_mat[lambda_ids[c],c]      # Deviance at optimal lambda for each cluster
  }
  lambda_vals = iClust_OS[[1]]$lambda[lambda_ids]

  # Order selected by deviance value (problem: maximum always at largest K index)
  # Alternative: set arbitrary threshold of 0.05: order selected when deviance doesn't increase by at least 0.05
  devs_shifted = c(NA,devs[-length(devs)])     # everything shifted to the right by 1 --> first elt = NA
  dev_inc = (devs-devs_shifted)/devs_shifted   # first elt = NA. 2nd elt: (dev2-dev1)/dev1, 3rd elt: (dev3-dev2)/dev2, ...

  # only one elt of K_search --> just choose the first K_search value
  if(all(is.na(dev_inc))){max_k=iCl_K_search[1]}

  if(any(dev_inc[-1]<0.05)){    # optimal cluster selected at index right before (less than 5% increase in POD) is observed
    max_k = iCl_K_search[which(dev_inc[-1]<0.05)[1]]
  } else{ #if no dev inc is <0.05, then last K selected
    max_k=iCl_K_search[length(iCl_K_search)]
  }
  max_lambda = lambda_vals[max_k]

  cat("Tuning complete.\n")
  cat(paste("Optimal K:",max_k,"+1\n"))

  fit <- iClusterPlus(dt1=t(cts),type="poisson",lambda=max_lambda,K=max_k,maxiter=10)
  end_iCl = as.numeric(Sys.time())

  cat(paste("Final fit complete. Total time elapsed: ",end_iCl-start_iCl,"s.\n",sep=""))

  cls_metrics = cluster.stats(d=NULL,fit$clusters,true_cls)

  #cat(cls_metrics)

  iCl_disc = if(is.matrix(fit$beta[[1]])){rowSums(fit$beta[[1]]==0)}else{fit$beta[[1]]==0}
  TPR = sum(iCl_disc & DEg_ID)/sum(DEg_ID)
  FPR = sum(iCl_disc & !DEg_ID)/sum(!DEg_ID)

  return(list(fit=fit,
              K=max_k+1,
              cls=fit$clusters,
              cls_metrics=cls_metrics,
              time_elap=(end_iCl-start_iCl),TPR=TPR,FPR=FPR))
}
run_HC=function(cts,idx,true_cls,K_search){
  library(FSCseq)
  start_HC = Sys.time()

  # Calculate distance metric to perform cluster validation, and for HC
  d = as.dist(1-cor(cts, method="spearman")) # distance metric on entire dataset

  cts=cts[idx,]  # prefilter
  set.seed(123)

  # gap statistic for OS
  results_hc=NbClust(data=t(log(cts+0.1)), diss=d, distance=NULL, min.nc=min(K_search), max.nc=max(K_search), method="average", index="gap")

  K = results_hc$Best.nc[1]
  cat(paste("HC K: ",K, "\n"))
  cat("Tuning complete.\n")

  cls = cutree(hclust(d,method="average"),K)
  end_HC = Sys.time()

  cls_metrics = cluster.stats(d=NULL,cls,true_cls,compareonly=T)

  return(list(fit=results_hc,
              K=K,
              cls=cls,
              cls_metrics=cls_metrics,
              time_elap=(end_HC-start_HC)))
}
run_KM=function(cts,true_cls,K_search){
  start_KM = Sys.time()
  if(K_search[1]==1){K_search=K_search[-1]}  # iCl can't compare K=1 --> remove this

  pam_fit = list()
  sil.vals = rep(NA,length(K_search))
  for(c in 1:length(K_search)){
    pam_fit[[c]] = pam(t(log(cts+0.1)),K_search[c])
    sil.vals[c] = pam_fit[[c]]$silinfo$avg.width
    cat(paste("KM order select avg silhouette: ", sil.vals[c],"\n"))
  }
  K = K_search[which.max(sil.vals)]
  cat(paste("KM",K, "\n"))
  cls = pam(t(log(cts+0.1)),K)$cluster

  end_KM = Sys.time()

  cls_metrics = cluster.stats(d=NULL,cls,true_cls,compareonly=T)

  return(list(fit=pam_fit[[which.max(sil.vals)]],
              K=K,
              cls=cls,
              cls_metrics=cls_metrics,
              time_elap=(end_KM-start_KM)))
}
run_NBMB=function(cts,true_cls,K_search){
  start_NBMB=Sys.time()
  if(K_search[1]==1){K_search=K_search[-1]}  # iCl can't compare K=1 --> remove this
  cts=round(cts,0)
  fit=NB.MClust(Count=t(cts),K=K_search)
  end_NBMB=Sys.time()

  cls_metrics = cluster.stats(d=NULL,fit$cluster,true_cls,compareonly=T)
  #cat(cls_metrics)

  return(list(fit=fit,
              K=fit$K,
              cls=fit$cluster,
              cls_metrics=cls_metrics,
              time_elap=(end_NBMB-start_NBMB)))
}
run_MCs=function(cts,idx,true_cls,K_search,dds){
  # output lMC, vMC, and rMC
  # Mclust (log, variance-stabilizing, and rlog transforms)

  start_lMC = Sys.time()
  fit=Mclust(t(log(cts+0.1)),G=K_search)
  cat("mclust log transform K: ",fit$G,"\n")
  end_lMC = Sys.time()
  cls_metrics = cluster.stats(d=NULL,fit$classification,true_cls,compareonly=T)
  #cat(cls_metrics)
  lMC=list(K=fit$G,
           cls=fit$classification,
           cls_metrics=cls_metrics,
           time_elap=(end_lMC-start_lMC))
  rm("fit")

  start_vMC = Sys.time()
  vsd_cts=NA
  try(vsd_cts <- assay(DESeq2::varianceStabilizingTransformation(dds))[idx,])
  # transform sometimes gives error in lfproc --> subset y (original cts) to idx genes in dds first (prevent low count/noisy genes)
  if(all(is.na(vsd_cts))){
    vsd_cts = assay(DESeq2::varianceStabilizingTransformation(dds,fitType="mean"))[idx,]
  }
  fit=Mclust(t(vsd_cts),G=K_search)
  rm("vsd_cts")
  cat("mclust vsd transform K: ",fit$G,"\n")
  end_vMC = Sys.time()
  cls_metrics = cluster.stats(d=NULL,fit$classification,true_cls,compareonly=T)
  #cat(cls_metrics)
  vMC=list(K=fit$G,
           cls=fit$classification,
           cls_metrics=cls_metrics,
           time_elap=(end_vMC-start_vMC))
  rm("fit")

  start_rMC = Sys.time()
  rld_cts=NA
  try(rld_cts <- assay(DESeq2::rlogTransformation(dds))[idx,])
  if(all(is.na(rld_cts))){
    rld_cts = assay(DESeq2::rlogTransformation(dds,fitType="mean"))[idx,]
  }
  fit=Mclust(t(rld_cts),G=K_search)
  rm("rld_cts")
  cat("mclust rld transform K: ",fit$G,"\n")
  end_rMC = Sys.time()
  cls_metrics = cluster.stats(d=NULL,fit$classification,true_cls,compareonly=T)
  rMC=list(K=fit$G,
           cls=fit$classification,
           cls_metrics=cls_metrics,
           time_elap=(end_rMC-start_rMC))
  rm("fit")

  return(list(lMC=lMC,vMC=vMC,rMC=rMC))
}
run_FSCseq_predict=function(sim.dat,FSC_summary){
  library(FSCseq)
  cts_pred=sim.dat$cts_pred; cls_pred=sim.dat$cls_pred
  fit=FSC_summary$fit; idx=FSC_summary$idx

  FSCseq_predict_results = FSCseq::FSCseq_predict_workflow(fit=fit,cts=sim.dat$cts,cts_pred=cts_pred,idx=idx)

  res=FSCseq_predict_results$results
  pARI = adjustedRandIndex(res$clusters,cls_pred)

  cls_metrics=cluster.stats(d=NULL,res$clusters,cls_pred,compareonly=T)
  return(list(fit=res,
              cls=res$clusters,
              pARI=pARI,
              cls_metrics=cls_metrics))
}
