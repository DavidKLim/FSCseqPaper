# Creates R scripts, and submits them as batch jobs in slurm

runScripts_FSC_init = function(ncores=10, dataset, med_filt=500, MAD_filt=50,
                      covariates='T',K_search=c(2:8)){
  lambda0=0.05; alpha0=0.01
  dname = sprintf("%s_med%d_MAD%d",dataset,med_filt,MAD_filt)
  ifelse(!dir.exists("./RealData/scripts"),
         dir.create("./RealData/scripts"),
         FALSE)
  ifelse(!dir.exists(sprintf("./RealData/scripts/%s",dname)),
         dir.create(sprintf("./RealData/scripts/%s",dname)),
         FALSE)

  for(c in 1:length(K_search)){
    cmd = rep(0, 3)
    cmd[1] = "unlink('.RData') \n setwd('./RealData') \n source('runRealAnalyses.R') \n"
    cmd[2] = sprintf("res = runRealAnalyses_FSC_init(ncores=%d,dataset='%s',med_filt=%d,MAD_filt=%d,
                                            K=%d,lambda=%f,alpha=%f,covariates=%s)\n",
                     ncores,dataset,med_filt,MAD_filt,
                     K_search[c],lambda0,alpha0,covariates)

    fname = sprintf("init%d_%f_%f_covars%s",
                    K_search[c],lambda0,alpha0,covariates)

    out = sprintf("./RealData/scripts/%s/%s",dname,fname)
    out2 = sprintf("./RealData/%s/%s.out",
                   dname,fname)

    if(file.exists(out2)){next}

    cmd[3] = sprintf("save(res, file = './%s/%s.out')", dname,fname)
    cmdf = paste(cmd, collapse = "")
    write.table(cmdf, file = out, col.names = F, row.names = F, quote = F)

    # memory usage
    mem = ncores*(7000*(dataset=="BRCA_full") + 5000*(dataset=="BRCA_pure"))
    # number of hours adjusted by dataset and whether to adjust for covariates
    hrs = (K_search[c]/ncores)*((40*(covariates=="F") + 200*(covariates=="T"))*(dataset=="BRCA_full") + (20*(covariates=="F") + 80*(covariates=="T"))*(dataset=="BRCA_pure"))

    run = sprintf("sbatch -p general -N 1 -n %d --mem=%d -t %d:00:00 -J %s_%s --wrap='R CMD BATCH %s'",
                  ncores,mem,round(hrs,0),dname,fname,out)
    print(out)
    Sys.sleep(1)
    system(run)
  }
}

# FSC_tune is designed to automatically save previous output --> can pick up where it left off if time runs out.
runScripts_FSC_tune = function(ncores=1,dataset,med_filt=500,MAD_filt=50,covariates='T',K_search=c(2:8),
                      lambda_search='seq(0.25,5,0.25)',alpha_search='c(0.01,seq(0.05,0.5,0.05))',hrs=24){
  ## RUN tune_order() and determine_order() FIRST, THEN USE THOSE RESULTS "<dataset>/init_res.out" TO INITIALIZE EM ##
  dname = sprintf("%s_med%d_MAD%d",dataset,med_filt,MAD_filt)
  ifelse(!dir.exists("./RealData/scripts"),
         dir.create("./RealData/scripts"),
         FALSE)
  ifelse(!dir.exists(sprintf("./RealData/scripts/%s",dname)),
         dir.create(sprintf("./RealData/scripts/%s",dname)),
         FALSE)

  for(c in 1:length(K_search)){
    if(K_search[c]==1){next}        # don't submit jobs for K=1
    cmd = rep(0, 2)
    cmd[1] = "unlink('.RData') \n setwd('./RealData') \n source('runRealAnalyses.R') \n"
    cmd[2] = sprintf("runRealAnalyses_FSC_tune(ncores=%d,dataset='%s',med_filt=%d,MAD_filt=%d,K=%d,
                           lambda_search=%s,alpha_search=%s,covariates=%s)\n",
                     ncores,dataset,med_filt,MAD_filt,K_search[c],lambda_search,alpha_search,covariates)

    ## edit from here
    fname = sprintf("joint%d_covars%s",
                    K_search[c],covariates)

    out = sprintf("./RealData/scripts/%s/%s",
                  dname,fname)

    cmdf = paste(cmd, collapse = "")
    write.table(cmdf, file = out, col.names = F, row.names = F, quote = F)

    # have to test mem/hrs on K=1-3 maybe, then extrapolate
    mem = ncores*(7000*(dataset=="BRCA_full") + 5000*(dataset=="BRCA_pure"))
    run = sprintf("sbatch -p general --mail-type='FAIL,END' --mail-user=deelim@live.unc.edu -N 1 -n %d --mem=%d -t %d:00:00 -J %s_%s --wrap='R CMD BATCH %s'",
                  ncores,mem,round(hrs,0),dname,fname,out)
    print(out)
    Sys.sleep(0.2)
    system(run)
  }
}

runScripts_MC=function(dataset,med_filt,MAD_filt,K_search='c(2:8)',hrs=24){
  dname = sprintf("%s_med%d_MAD%d",dataset,med_filt,MAD_filt)
  ifelse(!dir.exists("./RealData/scripts"),
         dir.create("./RealData/scripts"),
         FALSE)
  ifelse(!dir.exists(sprintf("./RealData/scripts/%s",dname)),
         dir.create(sprintf("./RealData/scripts/%s",dname)),
         FALSE)

  cmd = rep(0, 2)
  cmd[1] = "unlink('.RData') \n setwd('./RealData') \n source('runRealAnalyses.R') \n"
  cmd[2] = sprintf("runRealAnalyses_MC(dataset='%s',med_filt=%d,MAD_filt=%d,K_search=%s)\n",
                   dataset,med_filt,MAD_filt,K_search)

  ## edit from here
  fname="MCs"
  fname1 = sprintf("./RealData/%s/lMC_summary.out",dname)
  fname2 = sprintf("./RealData/%s/vMC_summary.out",dname)
  fname3 = sprintf("./RealData/%s/rMC_summary.out",dname)

  if(file.exists(fname1) & file.exists(fname2) & file.exists(fname3)){return("MC outputs exist. Jobs not submitted")}

  out = sprintf("./RealData/scripts/%s/%s",dname,fname)
  cmdf = paste(cmd, collapse = "")
  write.table(cmdf, file = out, col.names = F, row.names = F, quote = F)

  # have to test mem/hrs on K=1-3 maybe, then extrapolate
  mem = 30000    # MC requires lots of memory
  run = sprintf("sbatch -p general --mail-type='FAIL,END' --mail-user=deelim@live.unc.edu -N 1 -n 1 --mem=%d -t %d:00:00 -J %s_%s --wrap='R CMD BATCH %s'",
                mem,round(hrs,0),dname,fname,out)
  print(out)
  Sys.sleep(1)

  system(run)
}

runScripts_Others=function(dataset,med_filt,MAD_filt,K_search='c(2:8)',hrs=24){
  dname = sprintf("%s_med%d_MAD%d",dataset,med_filt,MAD_filt)
  ifelse(!dir.exists("./RealData/scripts"),
         dir.create("./RealData/scripts"),
         FALSE)
  ifelse(!dir.exists(sprintf("./RealData/scripts/%s",dname)),
         dir.create(sprintf("./RealData/scripts/%s",dname)),
         FALSE)
  cmd = rep(0, 2)
  cmd[1] = "unlink('.RData') \n setwd('./RealData') \n source('runRealAnalyses.R') \n"
  cmd[2] = sprintf("runRealAnalyses_others(dataset='%s',med_filt=%d,MAD_filt=%d,K_search=%s)\n",
                   dataset,med_filt,MAD_filt,K_search)

  ## edit from here
  fname="Others"
  fname1 = sprintf("./RealData/%s/HC_summary.out",dname)
  fname2 = sprintf("./RealData/%s/KM_summary.out",dname)
  fname3 = sprintf("./RealData/%s/NBMB_summary.out",dname)

  if(file.exists(fname1) & file.exists(fname2) & file.exists(fname3)){return("HC/KM/NBMB outputs exist. Jobs not submitted")}

  out = sprintf("./RealData/scripts/%s/%s",dname,fname)
  cmdf = paste(cmd, collapse = "")
  write.table(cmdf, file = out, col.names = F, row.names = F, quote = F)

  # have to test mem/hrs on K=1-3 maybe, then extrapolate
  mem = 8000
  run = sprintf("sbatch -p general --mail-type='FAIL,END' --mail-user=deelim@live.unc.edu -N 1 -n 1 --mem=%d -t %d:00:00 -J %s_%s --wrap='R CMD BATCH %s'",
                mem,round(hrs,0),dname,fname,out)
  print(out)
  Sys.sleep(1)

  system(run)
}

runScripts_iCl=function(ncores=25,dataset,med_filt,MAD_filt,K_search='c(2:8)',hrs=24){
  dname = sprintf("%s_med%d_MAD%d",dataset,med_filt,MAD_filt)
  ifelse(!dir.exists("./RealData/scripts"),
         dir.create("./RealData/scripts"),
         FALSE)
  ifelse(!dir.exists(sprintf("./RealData/scripts/%s",dname)),
         dir.create(sprintf("./RealData/scripts/%s",dname)),
         FALSE)
  cmd = rep(0, 2)
  cmd[1] = "unlink('.RData') \n setwd('./RealData') \n source('runRealAnalyses.R') \n"
  cmd[2] = sprintf("runRealAnalyses_iCl(ncores=%d,dataset='%s',med_filt=%d,MAD_filt=%d,K_search=%s)\n",
                   ncores,dataset,med_filt,MAD_filt,K_search)
  fname = "iCl"
  fname1 = sprintf("./RealData/%s/iCl_summary.out",dname)
  if(file.exists(fname1)){return("iCl output exists. Job not submitted")}

  out = sprintf("./RealData/scripts/%s/%s",dname,fname)

  cmdf = paste(cmd, collapse = "")
  write.table(cmdf, file = out, col.names = F, row.names = F, quote = F)


  mem = 10000

  run = sprintf("sbatch -p general --mail-type='END,FAIL' --mail-user=deelim@live.unc.edu -N 1 -n %d --mem=%d -t %d:00:00 -J %s_%s --wrap='R CMD BATCH %s'",
                ncores,mem,round(hrs,0),dname,fname,out)
  print(out)
  Sys.sleep(1)
  system(run)
}

