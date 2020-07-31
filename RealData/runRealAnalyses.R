runRealAnalyses_FSC_tune = function(ncores=1,dataset=c("BRCA_full","BRCA_pure"),med_filt=500,MAD_filt=50,
                                K,lambda_search,alpha_search,
                                covariates=T,nMB=1){
  library(FSCseq)
  dir_name="."
  load(sprintf("%s/%s_env.RData",dir_name,dataset))
  idx = (rowMeds>=med_filt) & (mads >= quantile(mads,MAD_filt/100))
  trace.folder=sprintf("%s/%s_med%d_MAD%d",dir_name,dataset,med_filt,MAD_filt)
  ifelse(!dir.exists(sprintf("%s",trace.folder)),
         dir.create(sprintf("%s",trace.folder)),
         FALSE)
  ifelse(!dir.exists(sprintf("%s/Diagnostics",trace.folder)),
         dir.create(sprintf("%s/Diagnostics",trace.folder)),
         FALSE)
  if(!covariates){X=NULL}    # user setting to not adjust for covariates in analysis


  load(sprintf("%s/init%d_0.050000_0.010000_covars%s.out",
               trace.folder,K,disp,method,substr(as.character(covariates),1,1)))
  res_K=res

  BICs=matrix(nrow=length(lambda_search)*length(alpha_search),ncol=3)
  time_elap=rep(NA,length(lambda_search)*length(alpha_search))
  index=1
  list_res=list()
  for(a in 1:length(alpha_search)){for(l in 1:length(lambda_search)){

    res.file=sprintf("%s/joint%d_%f_%f_covars%s.out",
                     trace.folder,K,lambda_search[l],alpha_search[a],substr(as.character(covariates),1,1))

    if(file.exists(res.file)){load(res.file)}else{
      if(l==1){
        # for new alpha (first value of lambda)
        res = FSCseq::FSCseq(ncores=ncores,X=X, y=cts[idx,], k=K, lambda=lambda_search[l], alpha=alpha_search[a],
                             size_factors=SF, norm_y=norm_y[idx,],true_clusters=as.numeric(as.factor(cls)), true_disc=NULL,
                             init_parms=T, init_cls=res_K$clusters, init_wts=res_K$wts, init_coefs=res_K$coefs, init_phi=res_K$phi,trace=T,
                             trace.file=sprintf("%s/Diagnostics/%d_%f_%f_covars%s.txt",
                                                trace.folder,K,lambda_search[l],alpha_search[a],substr(as.character(covariates),1,1)),
                             mb_size=floor(sum(idx)/nMB))
        # list_res[[c]]=res     # save the first EM run for each K
      } else if(l>1){
        # for subsequent value of lambda (load prev lambda val res)
        load(prev_res.file)
        res = FSCseq::FSCseq(ncores=ncores,X=X, y=cts[idx,], k=K, lambda=lambda_search[l], alpha=alpha_search[a],
                             size_factors=SF, norm_y=norm_y[idx,],true_clusters=as.numeric(as.factor(cls)), true_disc=NULL,
                             init_parms=T, init_cls=res$clusters, init_wts=res$wts, init_coefs=res$coefs, init_phi=res$phi,trace=T,
                             trace.file=sprintf("%s/Diagnostics/%d_%f_%f_covars%s.txt",
                                                trace.folder,K,lambda_search[l],alpha_search[a],substr(as.character(covariates),1,1)),
                             mb_size=floor(sum(idx)/nMB))
      }
    }

    # save res
    save(res,file=res.file)

    # save prev_res.file. loaded if l>1
    prev_res.file=res.file

    # save BICs/time_elap and print
    BICs[index,] = c(lambda_search[l],alpha_search[a],res$BIC)
    time_elap[index] = res$time_elap
    print(paste("(K,l,a)=(",K,",",BICs[index,1],",",BICs[index,2],"), BIC=",round(BICs[index,3],0),", time=",round(time_elap[index],1),"s.",
                sep=""))
    index=index+1
  }}

}


runRealAnalyses_FSC_init = function(ncores=1,dataset=c("BRCA_full","BRCA_pure"),med_filt=500,MAD_filt=50,
                           K, lambda=0.05, alpha=0.01, covariates=T, nMB=5){

  library(FSCseq)
  library(matrixStats)
  dir_name="."

  load(sprintf("%s/%s_env.RData",dir_name,dataset))
  # file name: "<dataset>_env.RData" contains real data pre-processed data/environment objects:
  ## cts (full counts g x n matrix)
  ## norm_y (normalized counts g x n matrix)
  ## SF (SFs length n vector)
  ## cls (length n, true cluster labels)
  ## X (design n x p matrix)
  ## mads (MAD values for each gene)
  ## rowMeds (median value for each gene)
  ## anno (list of annotations/colData for each sample (columns of cts). each elt of list is of length cts)

  idx = (rowMeds>=med_filt) & (mads >= quantile(mads,MAD_filt/100))
  trace.folder=sprintf("%s/%s_med%d_MAD%d",dir_name,dataset,med_filt,MAD_filt)
  ifelse(!dir.exists(sprintf("%s",trace.folder)),
         dir.create(sprintf("%s",trace.folder)),
         FALSE)
  ifelse(!dir.exists(sprintf("%s/Diagnostics",trace.folder)),
         dir.create(sprintf("%s/Diagnostics",trace.folder)),
         FALSE)

  if(!covariates){X=NULL}    # user setting to not adjust for covariates in analysis

  res = FSCseq::FSCseq(ncores=ncores, X=X, y=cts[idx,], k=K, lambda=lambda, alpha=alpha,
                       size_factors=SF, norm_y=norm_y[idx,],
                       true_clusters=as.numeric(as.factor(cls)), true_disc=NULL,
                       init_parms=F, init_cls=NULL, init_wts=NULL, init_coefs=NULL, init_phi=NULL,
                       n_rinits=1,method="CEM",trace=T,trace.file=sprintf("%s/Diagnostics/%d_%f_%f_covars%s.txt",
                                                                          trace.folder,K,lambda,alpha,substr(as.character(covariates),1,1)),
                       mb_size=floor(sum(idx)/nMB))

  return(res)
}

runRealAnalyses_iCl = function(ncores=1,dataset,med_filt=500,MAD_filt=50,
                               K_search=c(2:8)){
  library(iClusterPlus)
  dir_name="."

  load(sprintf("%s/%s_env.RData",dir_name,dataset))
  idx = (rowMeds>=med_filt) & (mads >= quantile(mads,MAD_filt/100))
  cts = round(norm_y[idx,],0)

  start_iCl = as.numeric(Sys.time())
  if(K_search[1]==1){K_search=K_search[-1]}  # iCl can't compare K=1 --> remove this
  iCl_K_search = K_search-1      # iCl_K = K - 1 --> K = iCl_K + 1

  cts=round(cts,0)       # coerce normalized counts to integer

  # Tuning
  iClust_OS <- list()
  for(c in 1:length(iCl_K_search)){
    cv.fit = tune.iClusterPlus(cpus=ncores,dt1=t(cts),type="poisson",K=iCl_K_search[c],n.lambda=25,maxiter=20)
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

  #cat(cls_metrics)

  iCl_summary=list(fit=fit,
                   K=max_k+1,
                   cls=fit$clusters,
                   time_elap=(end_iCl-start_iCl))

  iCl_file = sprintf("%s/%s_med%d_MAD%d/iCl_summary.out",dir_name,dataset,med_filt,MAD_filt)

  save(iCl_summary,file=iCl_file)

}


runRealAnalyses_MC = function(dataset,med_filt=500,MAD_filt=50,
                              K_search=c(2:8)){
  # output lMC, vMC, and rMC
  # Mclust (log, variance-stabilizing, and rlog transforms)
  library(mclust)
  dir_name="."

  load(sprintf("%s/%s_env.RData",dir_name,dataset))
  idx = (rowMeds>=med_filt) & (mads >= quantile(mads,MAD_filt/100))
  cts = round(norm_y[idx,],0)

  lMC_file = sprintf("%s/%s_med%d_MAD%d/lMC_summary.out",dir_name,dataset,med_filt,MAD_filt)
  if(!file.exists(lMC_file)){
    start_lMC = Sys.time()
    fit=Mclust(t(log(cts+0.1)),G=K_search)
    cat("mclust log transform K: ",fit$G,"\n")
    end_lMC = Sys.time()
    lMC_summary=list(fit=fit,K=fit$G,
                     cls=fit$classification,
                     time_elap=(end_lMC-start_lMC))
    save(lMC_summary,file=lMC_file)
    rm("fit")
    rm("lMC_summary")
    print(lMC_file)
  }

  library(FSCseq)
  load(sprintf("%s/%s_env.RData",dir_name,dataset))
  dds=processData(cts)$dds
  library(SummarizedExperiment)

  vMC_file = sprintf("%s/%s_med%d_MAD%d/vMC_summary.out",dir_name,dataset,med_filt,MAD_filt)
  if(!file.exists(vMC_file)){
    start_vMC = Sys.time()
    vsd_cts = assay(DESeq2::varianceStabilizingTransformation(dds))[idx,]
    fit=Mclust(t(vsd_cts),G=K_search)
    rm("vsd_cts")
    cat("mclust vsd transform K: ",fit$G,"\n")
    end_vMC = Sys.time()
    vMC_summary=list(fit=fit,K=fit$G,
                     cls=fit$classification,
                     time_elap=(end_vMC-start_vMC))
    save(vMC_summary,file=vMC_file)
    rm("fit")
    rm("vMC_summary")
    print(vMC_file)
  }

  rMC_file=sprintf("%s/%s_med%d_MAD%d/rMC_summary.out",dir_name,dataset,med_filt,MAD_filt)
  if(!file.exists(rMC_file)){
    start_rMC = Sys.time()
    rld_cts = assay(DESeq2::rlogTransformation(dds))[idx,]
    # save rld_cts b/c transformation takes very long time
    save(rld_cts,file=sprintf("%s/%s_med%d_MAD%d/rld_cts.out",dir_name,dataset,med_filt,MAD_filt))
    fit=Mclust(t(rld_cts),G=K_search)
    rm("rld_cts")
    cat("mclust rld transform K: ",fit$G,"\n")
    end_rMC = Sys.time()
    rMC_summary=list(fit=fit,K=fit$G,
                     cls=fit$classification,
                     time_elap=(end_rMC-start_rMC))
    save(rMC_summary,file=rMC_file)
    rm("fit")
    rm("rMC_summary")
    print(rMC_file)
  }

}

# HC, KM, NBMB
runRealAnalyses_others = function(dataset,med_filt=500,MAD_filt=50,
                                  K_search=c(2:8)){
  library(NbClust)
  library(NB.MClust)
  library(cluster)
  dir_name="."

  load(sprintf("%s/%s_env.RData",dir_name,dataset))

  idx = (rowMeds>=med_filt) & (mads >= quantile(mads,MAD_filt/100))
  cts = norm_y[idx,]

  d = as.dist(1-cor(cts, method="spearman"))

  # HC
  start_HC = Sys.time()
  set.seed(123)
  results_hc=NbClust(data=t(log(cts+0.1)), diss=d, distance=NULL, min.nc=min(K_search), max.nc=max(K_search), method="average", index="gap")
  K = results_hc$Best.nc[1]
  cat(paste("HC K: ",K, "\n"))
  cat("Tuning complete.\n")
  cls = cutree(hclust(d,method="average"),K)
  end_HC = Sys.time()
  HC_summary=list(fit=results_hc,
                  K=K,
                  cls=cls,
                  time_elap=(end_HC-start_HC))

  save(HC_summary,file=sprintf("%s/%s_med%d_MAD%d/HC_summary.out",dir_name,dataset,med_filt,MAD_filt))


  # KM
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
  KM_summary=list(fit=pam_fit[[which.max(sil.vals)]],
                  K=K,
                  cls=cls,
                  time_elap=(end_KM-start_KM))
  save(KM_summary,file=sprintf("%s/%s_med%d_MAD%d/KM_summary.out",dir_name,dataset,med_filt,MAD_filt))


  # NBMB
  start_NBMB=Sys.time()
  cts=round(cts,0)
  fit=NB.MClust(Count=t(cts),K=K_search)
  end_NBMB=Sys.time()
  NBMB_summary=list(fit=fit,
                    K=fit$K,
                    cls=fit$cluster,
                    time_elap=(end_NBMB-start_NBMB))
  save(NBMB_summary,file=sprintf("%s/%s_med%d_MAD%d/NBMB_summary.out",dir_name,dataset,med_filt,MAD_filt))

}

runRealAnalyses_NMF = function(dataset,med_filt=500,MAD_filt=50,
                                  K_search=c(2:8)){

  library(NMF)
  library(fpc)
  nmf.options(shared.memory=F)
  nmf.options(pbackend="seq")

  start_NMF=Sys.time()
  if(K_search[1]==1){K_search=K_search[-1]}  # iCl can't compare K=1 --> remove this

  cat("\nFitting default NMF model...\n")

  dir_name="."
  load(sprintf("%s/%s_env.RData",dir_name,dataset))
  idx = (rowMeds>=med_filt) & (mads >= quantile(mads,MAD_filt/100))
  cts = norm_y[idx,]

  all_fit = NMF::nmf(x=cts,rank=K_search)
  #all_fit1 = NMF::nmf(x=log(cts+0.1),rank=K_search)   # ??
  opt_k_index = which.max(all_fit$measures$silhouette.consensus)[1]
  opt_k = all_fit$measures$rank[opt_k_index]  # Highest average silhouette width, if >1 highest, lowest index selected.

  fit = all_fit$fit[[opt_k_index]]
  fit$clusters <- apply(coef(fit), 2, function(x)which.max(x)[1])

  cat("\nFit complete.\n")

  end_NMF=Sys.time()

  cls_metrics = cluster.stats(d=NULL,fit$clusters,true_cls,compareonly=T)
  #cat(cls_metrics)
  NMF_summary = list(fit=fit,
              K=opt_k,
              cls=fit$clusters,
              cls_metrics=cls_metrics,
              time_elap=(end_NMF-start_NMF))
  save(NMF_summary,file=sprintf("%s/%s_med%d_MAD%d/NMF_summary.out",dir_name,dataset,med_filt,MAD_filt))
  return(NMF_summary)

}
