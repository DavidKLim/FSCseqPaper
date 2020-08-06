library(FSCseq)

# varying params
K=c(2,4)
n_k=c(25,50)
LFCg=c(1,2)
pDEg=c(0.025,0.05)
beta0=c(8,12)
phi0=c(0.15,0.35,0.50)

simAllData = function(K,n_k,LFCg,pDEg,beta0,phi0,B,LFCb,nsims){
  # fixed params
  g=10000; pDEb=0.5
  sigma_g=0.1; sigma_b=0   # Added Gaussian (nonparametric) noise
  n_pred=25
  save_dir = sprintf("Simulations/%f_%f/B%d_LFCb%d/BatchPred",sigma_g,sigma_b,B,LFCb)
  for(a in 1:length(K)){for(b in 1:length(n_k)){for(c in 1:length(LFCg)){for(d in 1:length(pDEg)){for(e in 1:length(beta0)){for(f in 1:length(phi0)){
    # Simulates data into "./Simulations/<sigma_g>_<sigma_b>/<K>_<n>_<LFCg>_<pDEg>_<beta0>_<phi0>_sim<>_data.RData"
    FSCseq::simulateData(K=K[a], B=B, g=g, n=K[a]*n_k[b], LFCg=LFCg[c], pDEg=pDEg[d], sigma_g=sigma_g,
                 LFCb=LFCb, pDEb=pDEb, sigma_b=sigma_b, beta0=beta0[e], phi0=phi0[f],
                 nsims=nsims, n_pred=n_pred, sim_batch_pred=TRUE, save_dir=save_dir,save_file=T)
  }}}}}}
}

# # no batch
# B=1; LFCb=0
# simAllData(K=2,n_k,LFCg,pDEg,beta0,phi0,B,LFCb,25)
# simAllData(K=4,n_k,LFCg,pDEg,beta0=8,phi0,B,LFCb,25)
# simAllData(K=4,n_k,LFCg,pDEg,beta0=12,phi0,B,LFCb,100)    # main results (Table 1): simulated 100 datasets here

# batch with batch effect=2
B=2; LFCb=2
simAllData(K,n_k=50,LFCg=2,pDEg=0.05,beta0=12,phi0=0.35,B,LFCb,25)

# batch with batch effect=3
B=2; LFCb=3
simAllData(K,n_k=50,LFCg=2,pDEg=0.05,beta0=12,phi0=0.35,B,LFCb,25)
