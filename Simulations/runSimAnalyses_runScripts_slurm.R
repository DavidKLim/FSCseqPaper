# Submits the R script of each simulation case via each comparative method as a separate batch job on slurm.
# Can alternatively change 'run' to run interactively, but not recommended for full set of simulation conditions
# sourced by 'submitJobs_slurm.R', which calls these functions to submit R batch jobs into slurm.
# Scripts should have already been created by 'runSimAnalyses_createScripts.R'

runFSCs=function(ncores=1,sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                 K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                 sim_index=c(1:25)){

  if(B==1){LFCb=0}
  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)

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
runMCs=function(sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                sim_index=c(1:25)){

  if(B==1){LFCb=0}
  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  for(a in 1:length(K)){for(b in 1:length(n)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta)){for(f in 1:length(phi)){for(g in 1:length(sim_index)){
    fname = sprintf("%d_%d_%f_%f_%f_%f_sim%d_MC",K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],sim_index[g])
    res.file = sprintf("%s/%s.out",dir_name1,fname)
    out_script=sprintf("%s/%s",dir_name2,fname)
    if(!file.exists(res.file)){
      mem = 10
      mins = 20*(n[b]==50) + 40*(n[b]==100) + 80*(n[b]==200)
      run = sprintf("sbatch -p general -N 1 -n 1 --mem=%dg -t %d:00 -J %s --wrap='R CMD BATCH %s'",mem,mins,fname,out_script)
      print(run)
      Sys.sleep(0.2)
      system(run)
    }
  }}}}}}}
}
runiCls=function(ncores=5,sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                 K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                 sim_index=c(1:25)){

  if(B==1){LFCb=0}
  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  for(a in 1:length(K)){for(b in 1:length(n)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta)){for(f in 1:length(phi)){for(g in 1:length(sim_index)){

    fname = sprintf("%d_%d_%f_%f_%f_%f_sim%d_iCl",K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],sim_index[g])
    res.file = sprintf("%s/%s.out",dir_name1,fname)
    out_script=sprintf("%s/%s",dir_name2,fname)
    if(!file.exists(res.file)){
      mem = (1500 + 1000*(n[b]==200))*ncores
      hrs = round((40 + 40*(n[b]==200))/ncores,0)
      run = sprintf("sbatch -p general -N 1 -n %d --mem=%d -t %d:00:00 -J %s --wrap='R CMD BATCH %s'",ncores,mem,hrs,fname,out_script)
      print(run)
      Sys.sleep(0.2)
      system(run)
    }
  }}}}}}}
}
runOthers=function(sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                   K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                   sim_index=c(1:25)){

  if(B==1){LFCb=0}
  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  for(a in 1:length(K)){for(b in 1:length(n)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta)){for(f in 1:length(phi)){for(g in 1:length(sim_index)){
    fname = sprintf("%d_%d_%f_%f_%f_%f_sim%d",K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],sim_index[g])
    res.file1 = sprintf("%s/%s_HC.out",dir_name1,fname)
    res.file2 = sprintf("%s/%s_KM.out",dir_name1,fname)
    res.file3 = sprintf("%s/%s_NBMB.out",dir_name1,fname)

    out_script=sprintf("%s/%s_others",dir_name2,fname)
    if(!file.exists(res.file1) | !file.exists(res.file2) | !file.exists(res.file3)){
      mem = 3
      mins = 60 + 60*(n[b]==200)    # 10-15 minutes. 20 for n=100, 30 for n=200
      run = sprintf("sbatch -p general -N 1 -n 1 --mem=%dg -t %d:00 -J %s --wrap='R CMD BATCH %s'",mem,mins,fname,out_script)
      print(run)
      Sys.sleep(0.2)
      system(run)
    }
  }}}}}}}
}
runNMFs=function(sigma_g=0.1,sigma_b=0, B=1, LFCb=0,
                 K=c(2,4), n=c(100,200), LFCg=c(1,2), pDEg=c(0.025,0.05), beta=c(8,12), phi=c(0.15,0.35,0.5),
                 sim_index=c(1:25)){
  if(B==1){LFCb=0}
  dir_name1=sprintf("./Simulations/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  dir_name2=sprintf("./Simulations/scripts/%f_%f/B%d_LFCb%d",sigma_g,sigma_b,B,LFCb)
  for(a in 1:length(K)){for(b in 1:length(n)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta)){for(f in 1:length(phi)){for(g in 1:length(sim_index)){
    fname = sprintf("%d_%d_%f_%f_%f_%f_sim%d",K[a],n[b],LFCg[c],pDEg[d],beta[e],phi[f],sim_index[g])
    res.file = sprintf("%s/%s_NMF.out",dir_name1,fname)

    out_script=sprintf("%s/%s_NMF",dir_name2,fname)
    if(!file.exists(res.file)){
      mem = 2000
      mins = 180 + 60*(n[b]==100) + 300*(n[b]==200)    # NEED TO PROFILE FOR N=50/100/200
      run = sprintf("sbatch -p general -N 1 -n 1 --mem=%d -t %d:00 -J %s_NMF --wrap='R CMD BATCH %s'",mem,mins,fname,out_script)
      print(run)
      Sys.sleep(0.2)
      system(run)
    }
  }}}}}}}
}
