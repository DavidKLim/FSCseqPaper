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


## correlation of cluster-disc genes, and selecting smallest model among bottom 10% of BICs. done on longleaf (slightly different scripts and directories)
# # gene-gene correlation plots
# load("BRCA_pure3_2_env.RData")
# filt_idx = rowMeds>=500 & mads>=quantile(mads,0.5)
# norm_y=norm_y[filt_idx,]
#load("Real_Data/BRCA_pure3_2_med500_MAD50/joint5_1.750000_0.010000_gene_CEM_covarsF.out")
#load("Real_Data/BRCA_pure3_2_med500_MAD50/joint3_0.750000_0.250000_gene_CEM_covarsT.out")
# discriminatory=res$discriminatory
# disc_norm_y = norm_y[discriminatory,]
# cors = cor(t(disc_norm_y))
# png("Real_Data/BRCA_pure3_2_med500_MAD50/CEM_covarsF_corPlot.png",res=120,width=900,height=900); corrplot(cors,tl.pos = "td",tl.cex = 0.2,mar=c(0,0,5,0),tl.col = 'black',type="upper",order="hclust",diag=F,main="Correlations of BRCA cluster-discriminatory genes, Not adjusted for Plate"); dev.off()
#
# # finding smallest model in bottom 10% of BIC
# K_search=2:8; lambda_search=seq(0.25,5,0.25); alpha_search=c(0.01,seq(0.05,0.5,0.05))
# setwd("/pine/scr/d/e/deelim/out/Real_Data/BRCA_pure3_2_med500_MAD50")
# npoints = length(K_search)*length(lambda_search)*length(alpha_search)
# BICs=matrix(NA,nrow=npoints,ncol=5)
# index=1
# colnames(BICs) = c("K","l","a","BIC","ndisc")
# for(c in 1:length(K_search)){for(l in 1:length(lambda_search)){for(a in 1:length(alpha_search)){
#   BICs[index,1:3] = c(K_search[c], lambda_search[l],alpha_search[a])
#   file.name=sprintf("joint%d_%f_%f_gene_CEM_covarsF.out", K_search[c],lambda_search[l],alpha_search[a])
#   if(file.exists(file.name)){load(file.name); BICs[index,4:5] = c(res$BIC, sum(res$discriminatory))}else{print(file.name)}
#   index=index+1
# }}}
#
# BICs[order(BICs[,4])[1:100,],]
#
# BICs[order(BICs[,4]),] [which.min(BICs[order(BICs[,4])[1:floor(npoints/10)],][,5]),]    # lowest ndisc from bottom 10% of BIC
# # covarT : 3          1.00       0.30 9789277       1324
# # covarF : 2.00       0.25       0.50 8704580.49    1169.00
