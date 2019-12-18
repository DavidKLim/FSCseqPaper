######### BASE ON FINALIZE.R, AND Rnw DOCUMENTS. DELETE THIS COMMENT
######### CREATE: Supp Table 1-2 (stab1.out & stab2.out: FSCseq full res from K=2 & K=4)
######### CREATE: Supp Table 3-4 (stab3.out & stab4.out: K=2, comp methods, uncovered K* & average ARI)
######### CREATE: Supp Table 5-6 (stab5.out & stab6.out: K=4, comp methods, uncovered K* & average ARI)
######### CREATE: Supp Table 7 (df_time.out, rename stab7.out)

#############
## Figures ##
#############

# Figure 1: TPR vs FPR, K=4
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

dir_name="./Simulations/Results"
load(sprintf("%s/df2.out",dir_name))
load(sprintf("%s/df4.out",dir_name))

all_sims_df2 = df2$all_sims_df
df2=df2$df
all_sims_df4 = df4$all_sims_df
df4=df4$df

all_sims_df2$n=as.factor(all_sims_df2$n)
all_sims_df2$LFCg=as.factor(all_sims_df2$LFCg)
all_sims_df2$pDEg=as.factor(all_sims_df2$pDEg)
all_sims_df2$beta=as.factor(all_sims_df2$beta)
all_sims_df2$phi=as.factor(all_sims_df2$phi)
all_sims_df4$n=as.factor(all_sims_df4$n)
all_sims_df4$LFCg=as.factor(all_sims_df4$LFCg)
all_sims_df4$pDEg=as.factor(all_sims_df4$pDEg)
all_sims_df4$beta=as.factor(all_sims_df4$beta)
all_sims_df4$phi=as.factor(all_sims_df4$phi)


all_sims_df4$OA_FSC = as.factor(all_sims_df4$OA_FSC==1)
all_sims_df2$OA_FSC = as.factor(all_sims_df2$OA_FSC==1)


p_K4_n100_LFC1 = ggplot(all_sims_df4[all_sims_df4$LFCg==1 & all_sims_df4$n==100,],aes(x=FPR_FSC,y=TPR_FSC,color=pDEg,shape=phi))+geom_point(aes(size=OA_FSC))+
  ggtitle("n=100, Low LFC") + xlab("FPR") + ylab("TPR") +
  scale_shape_manual(name=expression(phi[0]),values=c(0,2,4)) +
  scale_color_discrete(name = expression(p[DE]), labels = c("0.025", "0.05")) + ylim(c(0,1))+
  scale_size_discrete(name="Correct Order")
p_K4_n100_LFC2 = ggplot(all_sims_df4[all_sims_df4$LFCg==2 & all_sims_df4$n==100,],aes(x=FPR_FSC,y=TPR_FSC,color=pDEg,shape=phi))+geom_point(aes(size=OA_FSC))+
  ggtitle("n=100, Moderate LFC") + xlab("FPR") + ylab("TPR") +
  scale_shape_manual(name=expression(phi[0]),values=c(0,2,4)) +
  scale_color_discrete(name = expression(p[DE]), labels = c("0.025", "0.05")) + ylim(c(0,1))+
  scale_size_discrete(name="Correct Order")
p_K4_n200_LFC1 = ggplot(all_sims_df4[all_sims_df4$LFCg==1 & all_sims_df4$n==200,],aes(x=FPR_FSC,y=TPR_FSC,color=pDEg,shape=phi))+geom_point(aes(size=OA_FSC))+
  ggtitle("n=200, Low LFC") + xlab("FPR") + ylab("TPR") +
  scale_shape_manual(name=expression(phi[0]),values=c(0,2,4)) +
  scale_color_discrete(name = expression(p[DE]), labels = c("0.025", "0.05")) + ylim(c(0,1))+
  scale_size_discrete(name="Correct Order")
p_K4_n200_LFC2 = ggplot(all_sims_df4[all_sims_df4$LFCg==2 & all_sims_df4$n==200,],aes(x=FPR_FSC,y=TPR_FSC,color=pDEg,shape=phi))+geom_point(aes(size=OA_FSC))+
  ggtitle("n=200, Moderate LFC") + xlab("FPR") + ylab("TPR") +
  scale_shape_manual(name=expression(phi[0]),values=c(0,2,4)) +
  scale_color_discrete(name = expression(p[DE]), labels = c("0.025", "0.05")) + ylim(c(0,1))+
  scale_size_discrete(name="Correct Order")

mylegend1<-g_legend(p_K4_n100_LFC1)
p1 = arrangeGrob(p_K4_n100_LFC1 + theme(legend.position="none"),
                 p_K4_n100_LFC2 + theme(legend.position="none"),
                 p_K4_n200_LFC1 + theme(legend.position="none"),
                 p_K4_n200_LFC2 + theme(legend.position="none"),
                 nrow=2)
main=textGrob("TPR vs. FPR in Discovery of Cluster-discriminatory Genes",gp=gpar(fontsize=20,font=1))
png(sprintf("%s/Fig1.png",dir_name),width=800,height=800)
grid.arrange(p1,right=mylegend1,ncol=1,heights=c(10, 1),top=main)
dev.off()


####################################################################################################################

# Figure 2: ARI, OA barplots stratified K=2/K=4. all methods
df2_subs = df2[df2$n %in% c(50) & df2$LFCg %in% c(1) & df2$pDEg==0.05,]
df4_subs = df4[df4$n %in% c(100) & df4$LFCg %in% c(1) & df4$pDEg==0.05,]

library(reshape2)
# violin plots (Fig3 revised)
ARIs2=df2_subs[,paste("ARI_",c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC"),sep="")]
OAs2=df2_subs[,paste("OA_",c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC"),sep="")]
ARIs4=df4_subs[,paste("ARI_",c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC"),sep="")]
OAs4=df4_subs[,paste("OA_",c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC"),sep="")]

df_ARIs2=melt(ARIs2)
df_ARIs2$variable=substring(df_ARIs2$variable,5)
df_OAs2=melt(OAs2)
df_OAs2$variable=substring(df_OAs2$variable,4)

df_ARIs4=melt(ARIs4)
df_ARIs4$variable=substring(df_ARIs4$variable,5)
df_OAs4=melt(OAs4)
df_OAs4$variable=substring(df_OAs4$variable,4)

df_violin2=rbind(df_OAs2,df_ARIs2)
df_violin2$metric=c(rep("OA",nrow(df_OAs2)),rep("ARI",nrow(df_ARIs2)))
df_violin4=rbind(df_OAs4,df_ARIs4)
df_violin4$metric=c(rep("OA",nrow(df_OAs4)),rep("ARI",nrow(df_ARIs4)))

colnames(df_violin4)[1]="method";  colnames(df_violin2)[1]="method"
df_violin2$metric=factor(df_violin2$metric,levels=c("OA","ARI"))
df_violin4$metric=factor(df_violin4$metric,levels=c("OA","ARI"))
df_violin2$method = factor(df_violin2$method,levels=c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC"))
df_violin4$method = factor(df_violin4$method,levels=c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC"))

p1=ggplot(df_violin2,aes(x=method,y=value,fill=metric,color=metric)) + ylab("")+
  geom_violin() + ggtitle(expression(K[true] ~ '=2')) + geom_point(pch = 21, position = position_jitterdodge(jitter.width=0.25,jitter.height=0.001))
p2=ggplot(df_violin4,aes(x=method,y=value,fill=metric,color=metric)) + ylab("")+
  geom_violin() + ggtitle(expression(K[true] ~ '=4')) + geom_point(pch = 21, position = position_jitterdodge(jitter.width=0.25,jitter.height=0.001))

mylegend<-g_legend(p1)
# png(sprintf("%s/Fig3_K2_LFC12_K4_LFC23.png",dir_name))
png(sprintf("%s/Fig2.png",dir_name),width=600,height=600)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"),
                         nrow=2),right=mylegend,ncol=1)
dev.off()

####################################################################################################################

############
## Tables ##
############

# Table 1: Main results (K=4, beta=12)
load(sprintf("%s/dfmain.out",dir_name))
dfmain=dfmain$df

metrics=c("K","n","LFCg","pDEg","beta","phi","K_FSC","ARI_FSC","ARI_pFSC","TPR_FSC","FPR_FSC")

dfmain_subs=dfmain[,metrics]

df=cbind(dfmain_subs[,c(2:4,6)],dfmain_subs[,c(7:11)])
colnames(df)=c("$\\mathbf{n}$","$\\mathbf{LFC}$","$\\mathbf{p_{DE}}$","$\\boldsymbol{\\phi}$","$K^*$","$ARI$","$pARI$","$TPR$","$FPR$")
save(df,file=sprintf("%s/Tab1.out",dir_name))

####################################################################################################################

# Table 2: Batch (n_k=50, LFCg=2, pDEg=0.05, beta=12, phi=0.35)
load(sprintf("%s/df2batch.out",dir_name)); load(sprintf("%s/df4batch.out",dir_name))
df2batch=df2batch$df
df4batch=df4batch$df

metrics_K = paste("K_",c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC"),sep="")
metrics_ARI = paste("ARI_",c("FSC","iCl","HC","KM","NBMB","lMC","vMC","rMC"),sep="")

K2s=unlist(df2batch[,metrics_K]); K4s=unlist(df4batch[,metrics_K])
ARI2s=unlist(df2batch[,metrics_ARI]); ARI4s=unlist(df4batch[,metrics_ARI])

methods=c("$FSC$","","iCl","",
          "HC","","KM","",
          "NBMB","","lMC","",
          "vMC","","rMC","")
df_batch=data.frame(methods,rep(c(2,3),times=length(metrics_K)/2),K2s,ARI2s,K4s,ARI4s)

colnames(df_batch)=c("","$\\gamma_0$","$K^*_2$","ARI","$K^*_4$","ARI")

save(df_batch,file=sprintf("%s/Tab2.out",dir_name))

####################################################################################################################

##########################
## Supplementary Tables ##
##########################

# Supp Tables 1 & 2
load(sprintf("%s/df2.out",dir_name));load(sprintf("%s/df4.out",dir_name))
df2=df2$df; df4=df4$df
df2=df2[df2$LFCg%in%c(1,2),]
df4=df4[df4$LFCg%in%c(1,2),]

metrics=c("n","LFCg","pDEg","beta","phi","K_FSC","OA_FSC","ARI_FSC","ARI_pFSC","TPR_FSC","FPR_FSC")

df2_FSC = df2[,metrics]
colnames(df2_FSC)=c("$\\mathbf{n}$","$\\mathbf{LFC}$","$\\mathbf{p_{DE}}$","$\\boldsymbol{\\beta}$","$\\boldsymbol{\\phi}$","$K^*$","OA","ARI","pARI","TPR","FPR")

df4_FSC = df4[,metrics]
colnames(df4_FSC)=c("$\\mathbf{n}$","$\\mathbf{LFC}$","$\\mathbf{p_{DE}}$","$\\boldsymbol{\\beta}$","$\\boldsymbol{\\phi}$","$K^*$","OA","ARI","pARI","TPR","FPR")

save(df2_FSC,file=sprintf("%s/sTab1.out",dir_name))
save(df4_FSC,file=sprintf("%s/sTab2.out",dir_name))
####################################################################################################################

# Supp Tables 3-6
metrics_K=c("n","LFCg","pDEg","beta","phi","K_FSC","K_iCl","K_HC","K_KM","K_NBMB","K_lMC","K_vMC","K_rMC")
metrics_ARI=c("n","LFCg","pDEg","beta","phi","ARI_FSC","ARI_iCl","ARI_HC","ARI_KM","ARI_NBMB","ARI_lMC","ARI_vMC","ARI_rMC")

apply_highlights = function(x){
  res=paste("\\textcolor{red}{",x,"}",sep="")
  return(res)
}

# (K=2) subset relevant features
df_compare_K2 = df2[,metrics_K]
colnames(df_compare_K2)=c("$\\mathbf{n}$","$\\mathbf{LFC}$","$\\mathbf{p_{DE}}$","$\\boldsymbol{\\beta}$","$\\boldsymbol{\\phi}$","$K^*_{FSC}$","$K^*_{iCl}$","$K^*_{HC}$","$K^*_{KM}$","$K^*_{NBMB}$","$K^*_{lMC}$","$K^*_{vMC}$","$K^*_{rMC}$")
df_compare_K2[,6:length(metrics_K)] = format(round(df_compare_K2[,6:length(metrics_K)],2),nsmall=2)
df_compare_ARI2 = df2[,metrics_ARI]
colnames(df_compare_ARI2)=c("$\\mathbf{n}$","$\\mathbf{LFC}$","$\\mathbf{p_{DE}}$","$\\boldsymbol{\\beta}$","$\\boldsymbol{\\phi}$","$ARI_{FSC}$","$ARI_{iCl}$","$ARI_{HC}$","$ARI_{KM}$","$ARI_{NBMB}$","$ARI_{lMC}$","$ARI_{vMC}$","$ARI_{rMC}$")
df_compare_ARI2[,6:length(metrics_ARI)] = format(round(df_compare_ARI2[,6:length(metrics_ARI)],2),nsmall=2)

# (K=2) apply highlights to best performers
for(i in 1:nrow(df2)){
  x=abs(as.numeric(df_compare_K2[i,6:ncol(df_compare_K2)])-2)
  id_K = which(x==min(x)) + 5
  x=as.numeric(df_compare_ARI2[i,6:ncol(df_compare_ARI2)])
  id_ARI = which(x==max(x)) + 5
  df_compare_K2[i,id_K] = apply_highlights(df_compare_K2[i,id_K])
  df_compare_ARI2[i,id_ARI] = apply_highlights(df_compare_ARI2[i,id_ARI])
}

# (K=4) subset relevant features
df_compare_K4 = df4[,metrics_K]
colnames(df_compare_K4)=c("$\\mathbf{n}$","$\\mathbf{LFC}$","$\\mathbf{p_{DE}}$","$\\boldsymbol{\\beta}$","$\\boldsymbol{\\phi}$","$K^*_{FSC}$","$K^*_{iCl}$","$K^*_{HC}$","$K^*_{KM}$","$K^*_{NBMB}$","$K^*_{lMC}$","$K^*_{vMC}$","$K^*_{rMC}$")
df_compare_K4[,6:length(metrics_K)] = format(round(df_compare_K4[,6:length(metrics_K)],2),nsmall=2)
df_compare_ARI4 = df4[,metrics_ARI]
colnames(df_compare_ARI4)=c("$\\mathbf{n}$","$\\mathbf{LFC}$","$\\mathbf{p_{DE}}$","$\\boldsymbol{\\beta}$","$\\boldsymbol{\\phi}$","$ARI_{FSC}$","$ARI_{iCl}$","$ARI_{HC}$","$ARI_{KM}$","$ARI_{NBMB}$","$ARI_{lMC}$","$ARI_{vMC}$","$ARI_{rMC}$")
df_compare_ARI4[,6:length(metrics_ARI)] = format(round(df_compare_ARI4[,6:length(metrics_ARI)],2),nsmall=2)

# (K=4) apply highlights to best performers
for(i in 1:nrow(df4)){
  x=abs(as.numeric(df_compare_K4[i,6:ncol(df_compare_K4)])-4)
  id_K = which(x==min(x)) + 5
  x=as.numeric(df_compare_ARI4[i,6:ncol(df_compare_ARI4)])
  id_ARI = which(x==max(x)) + 5
  df_compare_K4[i,id_K] = apply_highlights(df_compare_K4[i,id_K])
  df_compare_ARI4[i,id_ARI] = apply_highlights(df_compare_ARI4[i,id_ARI])
}

save(df_compare_K2,file=sprintf("%s/sTab3.out",dir_name))
save(df_compare_ARI2,file=sprintf("%s/sTab4.out",dir_name))
save(df_compare_K4,file=sprintf("%s/sTab5.out",dir_name))
save(df_compare_ARI4,file=sprintf("%s/sTab6.out",dir_name))

####################################################################################################################

# Supp Table 7
load(sprintf("%s/df2.out",dir_name)); load(sprintf("%s/df4.out",dir_name))
df2=df2$df; df4=df4$df
df2$time_FSC = df2$time_FSC/60; df2$time_iCl = df2$time_iCl/60
df4$time_FSC = df4$time_FSC/60; df4$time_iCl = df4$time_iCl/60

time_ids=which(substr(colnames(df2),1,4)=="time")    # same for df2 and df4
time_n50=colMeans(rbind(df2[df2$n==50,time_ids],df4[df4$n==50,time_ids]))
time_n100=colMeans(rbind(df2[df2$n==100,time_ids],df4[df4$n==100,time_ids]))
time_n200=colMeans(rbind(df2[df2$n==200,time_ids],df4[df4$n==200,time_ids]))

df_time = cbind(time_n50,time_n100,time_n200)

rownames(df_time) = c("$FSC$","iCl*","HC","KM","NBMB","lMC","vMC","rMC")
colnames(df_time) = c("n=50","n=100","n=200")

save(df_time,file=sprintf("%s/sTab7.out",dir_name))
