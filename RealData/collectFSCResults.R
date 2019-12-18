collect_FSC(dataset=c("BRCA_full","BRCA_pure"),med_filt=500,MAD_filt=50,covariates='T',
            K_search=c(2:8),lambda_search=seq(0.25,5,0.25),alpha_search=c(0.01,seq(0.05,0.5,0.05))){

  setwd("./RealData")
  dname = sprintf("./%s_med%d_MAD%d",dataset,med_filt,MAD_filt)

  BICs=matrix(NA,nrow=length(K_search)*length(lambda_search)*length(alpha_search),ncol=4)
  index=1
  for(c in 1:length(K_search)){for(l in 1:length(lambda_search)){for(a in 1:length(alpha_search)){
    fname = sprintf("%s/joint%d_%f_%f_covars%s.out",
                    dname,K_search[c],lambda_search[l],alpha_search[a],covariates)
    if(file.exists(fname)){
      print(fname)
      load(fname)
    }else{
      index=index+1
      next
    }
    BICs[index,]=c(K_search[c],lambda_search[l],alpha_search[a],res$BIC)
    print(BICs[index,])
    index=index+1
  }}}
  max_id = which.min(BICs[,4])
  if(length(max_id)>1){
    warning("More than 1 combination of lambda/alpha yielded same BIC. Choosing first one")
    max_id=max_id[1]
  }
  max_k=BICs[max_id,1]; max_lambda=BICs[max_id,2]; max_alpha=BICs[max_id,3]
  print(paste(", max (K,lambda,alpha) = (",max_k,",",max_lambda,",",max_alpha,")",sep=""))

  file.name=sprintf("%s/joint%d_%f_%f_covars%s.out",dname,max_k,max_lambda,max_alpha,covariates)
  print(paste("file: ",file.name))

  file.copy(from=file.name,to=sprintf("%s/FSC_covars%s_summary.out",dname,covariates),recursive=T)
  return(file.name)
}

fname_pureF = collect_FSC(dataset="BRCA_pure",K_search=c(2:8),covariates="F")
fname_pureT = collect_FSC(dataset="BRCA_pure",K_search=c(2:8),covariates="T")
fname_fullF = collect_FSC(dataset="BRCA_full",K_search=c(2:15),covariates="F")
fname_fullT = collect_FSC(dataset="BRCA_full",K_search=c(2:15),covariates="T")


