# Create separate scripts for each method and for each sim_index to be run in batch. Submitted to slurm by 'runSimAnalyses_runScripts_slurm.R'
# sourced by 'submitJobs_slurm.R', which calls these functions to create the R script files which will be submitted
# into slurm as batch jobs

createScripts_FSC=function(sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                 K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                 K_search='c(2:6)', lambda_search='seq(0.25,5,0.25)', alpha_search='c(0.01,seq(0.05,0.5,0.05))',
                 sim_index=c(1:25),ncores=1,nMB=5){

  if(B==1){LFCb=0}

  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
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
      cmd[1] = "unlink('.RData') \n source('./Simulations/runSimAnalyses.R') \n"
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

createScripts_MC=function(sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                K_search='c(2:6)',sim_index=c(1:25)){

  if(B==1){LFCb=0}

  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir.create(dir_name1, recursive=TRUE, showWarnings = F)
  dir.create(dir_name2, recursive=TRUE, showWarnings = F)

  for(a in 1:length(K)){for(b in 1:length(n)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta)){for(f in 1:length(phi)){for(g in 1:length(sim_index)){

    fname = sprintf("%d_%d_%f_%f_%f_%f_sim%d_MC",K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],sim_index[g])
    res.file = sprintf("%s/%s.out",dir_name1,fname)

    out_script=sprintf("%s/%s",dir_name2,fname)
    # run only if outfile does not exist.
    if(!file.exists(res.file)){

      cmd = rep(0, 2)
      cmd[1] = "unlink('.RData') \n source('./Simulations/runSimAnalyses.R') \n"
      cmd[2] = sprintf("runSimAnalyses_MC(K=%d,n=%d,LFCg=%f,pDEg=%f,beta=%f,phi=%f,
                          sigma_g=%f,sigma_b=%f,B=%d,LFCb=%d,K_search=%s,sim_index=%d)\n",
                       K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],
                       sigma_g,sigma_b,B,LFCb,K_search,sim_index[g])
      cmdf = paste(cmd, collapse = "")
      write.table(cmdf, file = out_script, col.names = F, row.names = F, quote = F)
    }
  }}}}}}}
}

createScripts_iCl=function(ncores=5,sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                 K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                 K_search='c(2:6)',sim_index=c(1:25)){

  if(B==1){LFCb=0}
  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir.create(dir_name1, recursive=TRUE, showWarnings = F)
  dir.create(dir_name2, recursive=TRUE, showWarnings = F)

  for(a in 1:length(K)){for(b in 1:length(n)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta)){for(f in 1:length(phi)){for(g in 1:length(sim_index)){

    fname = sprintf("%d_%d_%f_%f_%f_%f_sim%d_iCl",K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],sim_index[g])
    res.file = sprintf("%s/%s.out",dir_name1,fname)
    # run only if outfile does not exist

    out_script=sprintf("%s/%s",dir_name2,fname)

    if(!file.exists(res.file)){

      cmd = rep(0, 2)
      cmd[1] = "unlink('.RData') \n source('./Simulations/runSimAnalyses.R') \n"
      cmd[2] = sprintf("runSimAnalyses_iCl(ncores=%d,K=%d,n=%d,LFCg=%f,pDEg=%f,beta=%f,phi=%f,
                          sigma_g=%f,sigma_b=%f,B=%d,LFCb=%d,K_search=%s,sim_index=%d)\n",
                       ncores,K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],
                       sigma_g,sigma_b,B,LFCb,K_search,sim_index[g])
      cmdf = paste(cmd, collapse = "")
      write.table(cmdf, file = out_script, col.names = F, row.names = F, quote = F)
    }
  }}}}}}}
}

createScripts_Others=function(sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                   K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                   K_search='c(2:6)',sim_index=c(1:25)){
  # includes HC, KM, and NBMB

  if(B==1){LFCb=0}
  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir.create(dir_name1, recursive=TRUE, showWarnings = F)
  dir.create(dir_name2, recursive=TRUE, showWarnings = F)

  for(a in 1:length(K)){for(b in 1:length(n)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta)){for(f in 1:length(phi)){for(g in 1:length(sim_index)){

    fname = sprintf("%d_%d_%f_%f_%f_%f_sim%d",K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],sim_index[g])
    res.file1 = sprintf("%s/%s_HC.out",dir_name1,fname)
    res.file2 = sprintf("%s/%s_KM.out",dir_name1,fname)
    res.file3 = sprintf("%s/%s_NBMB.out",dir_name1,fname)

    out_script=sprintf("%s/%s_others",dir_name2,fname)

    # run only if HC/KM/NBMB res files do not exist
    if(!file.exists(res.file1) | !file.exists(res.file2) | !file.exists(res.file3)){

      cmd = rep(0, 2)
      cmd[1] = "unlink('.RData') \n source('./Simulations/runSimAnalyses.R') \n"
      cmd[2] = sprintf("runSimAnalyses_others(K=%d,n=%d,LFCg=%f,pDEg=%f,beta=%f,phi=%f,
                          sigma_g=%f,sigma_b=%f,B=%d,LFCb=%d,K_search=%s,sim_index=%d)\n",
                       K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],
                       sigma_g,sigma_b,B,LFCb,K_search,sim_index[g])
      cmdf = paste(cmd, collapse = "")
      write.table(cmdf, file = out_script, col.names = F, row.names = F, quote = F)
    }
  }}}}}}}
}

createScripts_NMF=function(sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                           K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                           K_search='c(2:6)',sim_index=c(1:25)){

  if(B==1){LFCb=0}
  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir.create(dir_name1, recursive=TRUE, showWarnings = F)
  dir.create(dir_name2, recursive=TRUE, showWarnings = F)

  for(a in 1:length(K)){for(b in 1:length(n)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta)){for(f in 1:length(phi)){for(g in 1:length(sim_index)){

    fname = sprintf("%d_%d_%f_%f_%f_%f_sim%d",K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],sim_index[g])
    res.file = sprintf("%s/%s_NMF.out",dir_name1,fname)

    out_script=sprintf("%s/%s_NMF",dir_name2,fname)

    # run only if HC/KM/NBMB res files do not exist
    if(!file.exists(res.file)){

      cmd = rep(0, 2)
      cmd[1] = "unlink('.RData') \n source('./Simulations/runSimAnalyses.R') \n"
      cmd[2] = sprintf("runSimAnalyses_NMF(K=%d,n=%d,LFCg=%f,pDEg=%f,beta=%f,phi=%f,
                          sigma_g=%f,sigma_b=%f,B=%d,LFCb=%d,K_search=%s,sim_index=%d)\n",
                       K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],
                       sigma_g,sigma_b,B,LFCb,K_search,sim_index[g])
      cmdf = paste(cmd, collapse = "")
      write.table(cmdf, file = out_script, col.names = F, row.names = F, quote = F)
    }
  }}}}}}}
}
