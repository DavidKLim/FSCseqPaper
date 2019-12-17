# source functions that perform each comparative method
source("./Simulations/runSimAnalyses_functions.R")

## Wrapper functions for each comparative method ##

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

  dir_name=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)

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

runSimAnalyses_iCl=function(ncores=5,K=2,n=100,LFCg=1,pDEg=0.05,beta=8,phi=0.15,
                            sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K_search=c(2:6),sim_index=1){
  library(fpc)
  library(mclust)
  library(cluster)
  library(NbClust)
  library(FSCseq)
  library(NB.MClust)
  library(iClusterPlus)
  library(pryr)
  library(SummarizedExperiment)

  # use file.name to specify data.file to load, and trace.folder to create to store diagnostics
  file.name = sprintf("./Simulations/%f_%f/B%d_LFCb%d/%d_%d_%f_%f_%f_%f_sim%d",sigma_g,sigma_b,B,LFCb,K,n,LFCg,pDEg,beta,phi,sim_index)
  data.file = sprintf("%s_data.RData",file.name)

  load(data.file)
  # Input simulated data
  cts = sim.dat$cts
  cls = sim.dat$cls
  batch = sim.dat$batch
  B=sim.dat$sim_params$B
  g=sim.dat$sim_params$g
  DEg_ID = sim.dat$DEg_ID

  # Process simulated data
  processed.dat = processData(y=cts,geoMeans=NULL,med_filt=TRUE,MAD_filt=TRUE,med_thresh=100,MAD_quant_thresh=50)
  SF = processed.dat$size_factors
  norm_y = processed.dat$norm_y
  idx = processed.dat$idx       # length g, TRUE/FALSE

  rm("sim.dat"); rm("processed.dat")

  # mem_used() # 643 MB

  if(K_search[1]==1){
    K_search2=K_search[-1]      # for methods that can't search K=1 (iCl, KM, NBMB), omit
  }else{K_search2=K_search}

  res2=sprintf("%s_iCl.out",file.name)
  if(!file.exists(res2)){
    cat("\nRunning iCl analysis...\n")
    iCl_summary = run_iCl(cts=norm_y[idx,],true_cls=cls,K_search=K_search2,
                          ncores=ncores,n.lambda=25,DEg_ID=DEg_ID[idx])      ######### CHANGE THIS TO 25 #####
    save(iCl_summary,file=res2)
  } else{
    cat("\niCl analysis already run...\n")
  }
}

# MC: lMC, vMC, rMC
runSimAnalyses_MC=function(K=2,n=100,LFCg=1,pDEg=0.05,beta=8,phi=0.15,
                           sigma_g=0.1,sigma_b=0,B=1,LFCb=0,
                           K_search=c(2:6),sim_index=1){
  library(fpc)
  library(mclust)
  library(cluster)
  library(pryr)
  library(SummarizedExperiment)
  library(FSCseq)

  # use file.name to specify data.file to load, and trace.folder to create to store diagnostics
  file.name = sprintf("./Simulations/%f_%f/B%d_LFCb%d/%d_%d_%f_%f_%f_%f_sim%d",sigma_g,sigma_b,B,LFCb,K,n,LFCg,pDEg,beta,phi,sim_index)
  data.file = sprintf("%s_data.RData",file.name)

  # if sim has already been run (results exist in Simulations folder) --> don't run
  res.file = sprintf("%s_MC.out",file.name)
  if(file.exists(res.file)){
    print("Simulation MC has already been run")
    load(res.file)
    return(MC_summary)
  }

  load(data.file)
  # Input simulated data
  cts = sim.dat$cts
  cls = sim.dat$cls
  batch = sim.dat$batch
  B=sim.dat$sim_params$B
  g=sim.dat$sim_params$g
  DEg_ID = sim.dat$DEg_ID

  # Process simulated data
  processed.dat = processData(y=cts,geoMeans=NULL,med_filt=TRUE,MAD_filt=TRUE,med_thresh=100,MAD_quant_thresh=50)
  SF = processed.dat$size_factors
  norm_y = processed.dat$norm_y
  idx = processed.dat$idx       # length g, TRUE/FALSE
  dds=processed.dat$dds

  rm("sim.dat")
  rm("processed.dat")

  #### CLUSTERING ON SIMULATED DATASET ####
  cat("\nRunning transformed MC analyses (takes most memory)...\n")
  MC_summary = run_MCs(cts=norm_y[idx,],idx=idx,true_cls=cls,K_search=K_search,
                       dds=dds)

  save(MC_summary,file=res.file)
}

# "others" = HC, KM, NBMB
runSimAnalyses_others=function(K=2,n=100,LFCg=1,pDEg=0.05,beta=8,phi=0.15,
                               sigma_g=0.1,sigma_b=0,B=1,LFCb=0,K_search=c(2:6),sim_index=1){
  library(fpc)
  library(mclust)
  library(cluster)
  library(NbClust)
  library(FSCseq)
  library(NB.MClust)
  library(pryr)

  # use file.name to specify data.file to load, and trace.folder to create to store diagnostics
  file.name = sprintf("./Simulations/%f_%f/B%d_LFCb%d/%d_%d_%f_%f_%f_%f_sim%d",sigma_g,sigma_b,B,LFCb,K,n,LFCg,pDEg,beta,phi,sim_index)
  data.file = sprintf("%s_data.RData",file.name)

  load(data.file)
  # Input simulated data
  cts = sim.dat$cts
  cls = sim.dat$cls
  batch = sim.dat$batch
  B=sim.dat$sim_params$B
  g=sim.dat$sim_params$g
  DEg_ID = sim.dat$DEg_ID

  # Process simulated data
  processed.dat = processData(y=cts,geoMeans=NULL,med_filt=TRUE,MAD_filt=TRUE,med_thresh=100,MAD_quant_thresh=50)
  SF = processed.dat$size_factors
  norm_y = processed.dat$norm_y
  idx = processed.dat$idx       # length g, TRUE/FALSE

  rm("sim.dat")
  rm("processed.dat")

  # mem_used() # 643 MB

  if(K_search[1]==1){
    K_search2=K_search[-1]      # for methods that can't search K=1 (iCl, KM, NBMB), omit
  }else{K_search2=K_search}

  res4=sprintf("%s_KM.out",file.name)
  if(!file.exists(res4)){
    cat("\nRunning KM analysis...\n")
    KM_summary = run_KM(cts=norm_y[idx,],true_cls=cls,K_search=K_search2)
    save(KM_summary,file=res4)
  } else{
    cat("\nKM analysis already run...\n")
  }

  res5=sprintf("%s_NBMB.out",file.name)
  if(!file.exists(res5)){
    cat("\nRunning NBMB analysis...\n")
    NBMB_summary = run_NBMB(cts=norm_y[idx,],true_cls=cls,K_search=K_search2)
    save(NBMB_summary,file=res5)
  } else{
    cat("\nNBMB analysis already run...\n")
  }

  res3=sprintf("%s_HC.out",file.name)
  if(!file.exists(res3)){
    cat("\nRunning HC analysis...\n")
    HC_summary = run_HC(cts=norm_y,idx=idx,true_cls=cls,K_search=K_search)
    save(HC_summary,file=res3)
  } else{
    cat("\nHC analysis already run...\n")
  }
}
