library(FSCseq)
runCV=function(dataset,med_filt,MAD_filt,covariates="F"){
  load(sprintf("%s_env.RData",dataset))
  load(sprintf("./%s_med%d_MAD%d/FSC_covars%s_summary.out",
               dataset,med_filt,MAD_filt,covariates))

  filt_ids = (rowMeds>=med_filt) & (mads >= quantile(mads,(MAD_filt)/100))

  y=cts[filt_ids,]
  norm_y=norm_y[filt_ids,]

  disc_ids = res$discriminatory

  lambda=0
  alpha=0
  k=res$k

  ids_disc = res$discriminatory
  pred_cls = rep(NA,ncol(y))
  fits_train=list()
  fits_test=list()
  for(i in 1:ncol(y)){
    if(covariates=="F"){Xo=NULL}else{Xo=X[-i,]}
    fits_train[[i]] = FSCseq::FSCseq(X=Xo,y=y[ids_disc,-i], k=k,
                                     lambda=0,alpha=0,
                                     size_factors=SF[-i],
                                     norm_y=norm_y[ids_disc,-i],
                                     true_clusters=res$clusters[-i],
                                     init_cls=res$clusters[-i], init_wts=res$wts[,-i], n_rinits=0, maxit_inits=1,
                                     maxit_EM=1,method="EM",mb_size=sum(ids_disc))           # just M step, and then output

    test_SF = SF[i]
    test_true_cl = res$clusters[i]

    if(covariates=="F"){Xp=NULL}else{Xp=matrix(X[i,],nrow=1)}
    fits_test[[i]] = FSCseq::FSCseq_predict(X=Xp,fit=fits_train[[i]],cts_pred=matrix(y[ids_disc,i],ncol=1),SF_pred=test_SF)

    pred_cls[i]=fits_test[[i]]$clusters
    print(paste("Leave-one-out CV Sample",i,"Prediction&truth agreed:",pred_cls[i]==test_true_cl))
  }
  pred_acc=mean(pred_cls==res$clusters)

  save(list=c("pred_cls","fits_train","fits_test","pred_acc"),
       file=sprintf("./%s_med%d_MAD%d/FSC_covar%s_CV.RData",
                    dataset,med_filt,MAD_filt,covariates))

}
setwd("./RealData")
runCV(dataset="BRCA_pure",med_filt=500,MAD_filt=50,covariates="F")
