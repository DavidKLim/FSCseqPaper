
collectSims=function(sigma_g=0.1,sigma_b=0,B=1,LFCb=0,
                     K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                     compared_methods=c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC","pFSC"),
                     sim_index=c(1:25)){

  dir_name=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
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
      result=(summary$cls_metrics$cluster.number)
      if(is.null(result)){result=1}
    } else if(metric=="OA"){
      if(is.null(summary$cls_metrics$cluster.number)){result=0} else{
        result=((summary$cls_metrics$cluster.number==true.K)^2)
      }
    } else if(metric=="time"){
      result=(summary$time_elap)
    } else if(metric=="ARI"){
      result=(summary$cls_metrics$corrected.rand)
    }

    # else if(metric=="ASW"){
    #   result=(summary$cls_metrics$avg.silwidth)
    # } else if(metric=="VI"){
    #   result=(summary$cls_metrics$vi)
    # } else if(metric=="PG"){
    #   result=(summary$cls_metrics$pearsongamma)
    # } else if(metric=="dunn"){
    #   result=(summary$cls_metrics$dunn)
    # } else if(metric=="dunn2"){
    #   result=(summary$cls_metrics$dunn2)
    # } else if(metric=="entropy"){
    #   result=(summary$cls_metrics$entropy)
    # } else if(metric=="ch"){
    #   result=(summary$cls_metrics$ch)
    # } else if(metric=="bARI"){
    #   if(is.null(batch)){result=NA} else{result=adjustedRandIndex(summary$cls,batch)}
    # }

    if(is.null(result)){return(NA)} else{return(result)}
  }

  # index for each combination of sim params
  index=1
  all_sims_df = data.frame(matrix(nrow=0,ncol=ncol(df)))

  # ARI_sigmas_df = data.frame(matrix(nrow=0,ncol=length(compared_methods)))
  # OA_sigmas_df = data.frame(matrix(nrow=0,ncol=length(compared_methods[compared_methods!="pFSC"])))
  # OA_25_df=OA_sigmas_df; OA_75_df=OA_sigmas_df
  # K_sigmas_df = data.frame(matrix(nrow=0,ncol=length(compared_methods[compared_methods!="pFSC"])))
  # K_25_df=K_sigmas_df; K_75_df=K_sigmas_df
  # TPR_sigma=rep(NA,ncol(df)); TPR_25=rep(NA,ncol(df)); TPR_75=rep(NA,ncol(df))
  # FPR_sigma=rep(NA,ncol(df)); FPR_25=rep(NA,ncol(df)); FPR_75=rep(NA,ncol(df))

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

#     ARIs = matrix(unlist(all_sims_mat[,paste("ARI_",compared_methods,sep="")]),ncol=length(compared_methods))
#     # print(ARIs)
#     sds = apply(ARIs,2,sd)
#     ARI_sigmas_df=rbind(ARI_sigmas_df,sds)
#
#     OAs = matrix(unlist(all_sims_mat[,paste("OA_",compared_methods[compared_methods!="pFSC"],sep="")]),ncol=length(compared_methods[compared_methods!="pFSC"]))
#     Ks = matrix(unlist(all_sims_mat[,paste("K_",compared_methods[compared_methods!="pFSC"],sep="")]),ncol=length(compared_methods[compared_methods!="pFSC"]))
#
#     # print(OAs)
#     sds = apply(OAs,2,sd)
#     OAs25 = apply(OAs,2,function(x)quantile(x,0.25,na.rm=T))
#     OAs75 = apply(OAs,2,function(x)quantile(x,0.75,na.rm=T))
#
#     OA_sigmas_df=rbind(OA_sigmas_df,sds)
#     OA_25_df=rbind(OA_25_df,OAs25)
#     OA_75_df=rbind(OA_75_df,OAs75)
#
#     sds = apply(Ks,2,sd)
#     Ks25 = apply(Ks,2,function(x)quantile(x,0.25,na.rm=T))
#     Ks75 = apply(Ks,2,function(x)quantile(x,0.75,na.rm=T))
#     K_sigmas_df=rbind(K_sigmas_df,sds)
#     K_25_df=rbind(K_25_df,Ks25)
#     K_75_df=rbind(K_75_df,Ks75)
#
#     TPRs=all_sims_mat[,"TPR_FSC"]
#     FPRs=all_sims_mat[,"FPR_FSC"]
#
#     TPR_sigma[index]=sd(TPRs); FPR_sigma[index]=sd(FPRs)
#     TPR_25[index]=quantile(TPRs,0.25,na.rm=T); TPR_75[index]=quantile(TPRs,0.75,na.rm=T)
#     FPR_25[index]=quantile(FPRs,0.25,na.rm=T); FPR_75[index]=quantile(FPRs,0.75,na.rm=T)

    # average out every variable
    index = index+1
  }}}}}}

  # colnames(all_sims_df) = colnames(df)
  # colnames(ARI_sigmas_df)=paste("ARI_",compared_methods,sep="")
  # colnames(OA_sigmas_df)=paste("OA_",compared_methods[compared_methods!="pFSC"],sep="")
  # colnames(OA_25_df)=paste("OA_",compared_methods[compared_methods!="pFSC"],sep="")
  # colnames(OA_75_df)=paste("OA_",compared_methods[compared_methods!="pFSC"],sep="")
  # colnames(K_sigmas_df)=paste("K_",compared_methods[compared_methods!="pFSC"],sep="")
  # colnames(K_25_df)=paste("K_",compared_methods[compared_methods!="pFSC"],sep="")
  # colnames(K_75_df)=paste("K_",compared_methods[compared_methods!="pFSC"],sep="")
  # ARI_sigmas_df = cbind(df[,1:6],ARI_sigmas_df)
  # OA_sigmas_df = cbind(df[,1:6],OA_sigmas_df)
  # OA_25_df = cbind(df[,1:6],OA_25_df)
  # OA_75_df = cbind(df[,1:6],OA_75_df)
  # K_sigmas_df = cbind(df[,1:6],K_sigmas_df)
  # K_25_df = cbind(df[,1:6],K_25_df)
  # K_75_df = cbind(df[,1:6],K_75_df)


  return(list(df=df,all_sims_df=all_sims_df
              # ,ARI_sigmas_df=ARI_sigmas_df,
              # OA_sigmas_df=OA_sigmas_df,
              # OA_25_df=OA_25_df,OA_75_df=OA_75_df,
              # K_sigmas_df=K_sigmas_df,
              # K_25_df=K_25_df,K_75_df=K_75_df,
              # TPR_sigma=TPR_sigma, FPR_sigma=FPR_sigma,
              # TPR_25=TPR_25,TPR_75=TPR_75, FPR_25=FPR_25,FPR_75=FPR_75
              )) # return result to evaluate interactively
}

res.dir="./Simulations/Results"
ifelse(!dir.exists(res.dir),dir.create(res.dir),FALSE)

# collect main results
df2=collectSims(sigma_g=0.1,sigma_b=0,B=1,LFCb=0,
                K=c(2), n=c(50,100), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                compared_methods=c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC","pFSC"),
                sim_index=c(1:25))
save(df2,file=sprintf("%s/df2.out",res.dir))

df4=collectSims(sigma_g=0.1,sigma_b=0,B=1,LFCb=0,
                K=c(4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                compared_methods=c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC","pFSC"),
                sim_index=c(1:25))
save(df4,file=sprintf("%s/df4.out",res.dir))

dfmain=collectSims(sigma_g=0.1,sigma_b=0,B=1,LFCb=0,
                K=c(4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(12), phi=c(0.15,0.35,0.5),
                compared_methods=c("FSC","pFSC"),
                sim_index=c(1:100))
save(dfmain,file=sprintf("%s/dfmain.out",res.dir))


# collect batch results
df2batch2=collectSims(sigma_g=0.1,sigma_b=0,B=2,LFCb=2,
                     K=c(2), n=c(100), LFCg=c(2), pDEg=c(0.05), beta=c(12), phi=c(0.35),
                     compared_methods=c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC","pFSC"),
                     sim_index=c(1:25))
df2batch2$df$gamma=2;df2batch2$all_sims_df$gamma=2
df4batch2=collectSims(sigma_g=0.1,sigma_b=0,B=2,LFCb=2,
                     K=c(4), n=c(200), LFCg=c(2), pDEg=c(0.05), beta=c(12), phi=c(0.35),
                     compared_methods=c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC","pFSC"),
                     sim_index=c(1:25))
df4batch2$df$gamma=2;df4batch2$all_sims_df$gamma=2
df2batch=list(df=rbind(df2batch2$df,df2batch3$df),all_sims_df=rbind(df2batch2$all_sims_df,df2batch3$all_sims_df))

df2batch3=collectSims(sigma_g=0.1,sigma_b=0,B=2,LFCb=3,
                      K=c(2), n=c(100), LFCg=c(2), pDEg=c(0.05), beta=c(12), phi=c(0.35),
                      compared_methods=c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC","pFSC"),
                      sim_index=c(1:25))
df2batch3$df$gamma=3;df2batch3$all_sims_df$gamma=3
df4batch3=collectSims(sigma_g=0.1,sigma_b=0,B=2,LFCb=3,
                      K=c(4), n=c(200), LFCg=c(2), pDEg=c(0.05), beta=c(12), phi=c(0.35),
                      compared_methods=c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC","pFSC"),
                      sim_index=c(1:25))
df4batch3$df$gamma=3;df4batch3$all_sims_df$gamma=3
df4batch=list(df=rbind(df4batch2$df,df4batch3$df),all_sims_df=rbind(df4batch2$all_sims_df,df4batch3$all_sims_df))

save(df2batch,file=sprintf("%s/df2batch.out",res.dir)); save(df4batch,file=sprintf("%s/df4batch.out",res.dir))



