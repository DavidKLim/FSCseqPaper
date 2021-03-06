# Replicate Tables/Figures pertaining to BRCA_pure and BRCA_full from FSCseq paper
setwd("./RealData")
tabulateResults = function(dataset=c("BRCA_pure","BRCA_full"),med_filt=500,MAD_filt=50){
  library(fpc)
  res_dir=sprintf("./%s_med%d_MAD%d",dataset,med_filt,MAD_filt)

  FSC.file_F=sprintf("%s/FSC_covarsF_summary.out",res_dir)
  FSC.file_T=sprintf("%s/FSC_covarsT_summary.out",res_dir)

  FSC2.file_F=sprintf("%s/FSC2_covarsF_summary.out",res_dir)   ###
  FSC2.file_T=sprintf("%s/FSC2_covarsT_summary.out",res_dir)   ###

  iCl.file=sprintf("%s/iCl_summary.out",res_dir)
  HC.file=sprintf("%s/HC_summary.out",res_dir)
  KM.file=sprintf("%s/KM_summary.out",res_dir)
  NMF.file=sprintf("%s/NMF_summary.out",res_dir)
  NBMB.file=sprintf("%s/NBMB_summary.out",res_dir)
  lMC.file=sprintf("%s/lMC_summary.out",res_dir)
  vMC.file=sprintf("%s/vMC_summary.out",res_dir)
  rMC.file=sprintf("%s/rMC_summary.out",res_dir)

  load(FSC.file_F); resF=res; load(FSC.file_T); resT=res

  load(FSC2.file_F); resF2=res; load(FSC2.file_T); resT2=res   ###

  load(iCl.file);load(HC.file);load(KM.file);load(NBMB.file);load(NMF.file);load(lMC.file);load(vMC.file);load(rMC.file)
  load(sprintf("%s_env.RData",dataset))

  cls=as.numeric(as.factor(cls))
  d = as.dist(1-cor(norm_y, method="spearman"))
  anno_summary=list(cls=cls,K=length(unique(cls)))
  FSC_summary=list(cls=resF$clusters,K=length(unique(resF$clusters)))
  FSCadj_summary=list(cls=resT$clusters,K=length(unique(resT$clusters)))

  FSC10_summary=list(cls=resF2$clusters,K=length(unique(resF2$clusters)))   ###
  FSCadj10_summary=list(cls=resT2$clusters,K=length(unique(resT2$clusters)))   ###

  # summaries = c("anno_summary","FSC_summary","FSCadj_summary","iCl_summary",
  #               "HC_summary","KM_summary","NBMB_summary","NMF_summary","lMC_summary","vMC_summary","rMC_summary")
  summaries = c("anno_summary","FSC_summary","FSC10_summary","FSCadj_summary","FSCadj10_summary","iCl_summary",
                "HC_summary","KM_summary","NBMB_summary","NMF_summary","lMC_summary","vMC_summary","rMC_summary")    ###

  metrics=c("K","RI","NID","NVI","NMI","ARI")
  #####RI(up),NID(down),NVI(down),NMI(up),ARI(up),ASW(up),VI(down),PG(up),dunn(up),dunn2(up),entropy(down),ch(up),wb.ratio(down),sindex(up)
  library(aricode)
  fill.row=function(d,summary){
    metrics=c("K","RI","NID","NVI","NMI","ARI")
    summary_obj = eval(parse(text=summary))
    res = cluster.stats(d,summary_obj$cls,cls)
    res2 = clustComp(summary_obj$cls,cls)
    res_vec = c(res$cluster.number,res2$RI,res2$NID,res2$NVI,res2$NMI,res2$ARI)
    names(res_vec) = metrics
    return(res_vec)
  }

  df=data.frame(matrix(nrow=length(summaries),ncol=length(metrics)))
  colnames(df)=metrics
  rownames(df)=sapply(strsplit(summaries,"_"),function(x) x[1])

  for(i in 1:length(summaries)){
    df[i,] = fill.row(d,summaries[i])
  }

  if(dataset == "BRCA_full"){
    resI=cluster.stats(d,as.numeric(as.factor(anno$cls_sigCL_intrinsic)),cls)
    resI2=clustComp(as.numeric(as.factor(anno$cls_sigCL_intrinsic)),cls)
    resU=cluster.stats(d,as.numeric(as.factor(anno$cls_sigCL_unsupervised)),cls)
    resU2=clustComp(as.numeric(as.factor(anno$cls_sigCL_unsupervised)),cls)
    res_vecI = c(resI$cluster.number,resI2$RI,resI2$NID,resI2$NVI,resI2$NMI,resI2$ARI,resI$avg.silwidth,resI$vi,resI$pearsongamma,resI$dunn,resI$dunn2,
                 resI$entropy,resI$ch,resI$wb.ratio,resI$sindex)
    res_vecU = c(resU$cluster.number,resU2$RI,resU2$NID,resU2$NVI,resU2$NMI,resU2$ARI,resU$avg.silwidth,resU$vi,resU$pearsongamma,resU$dunn,resU$dunn2,
                 resU$entropy,resU$ch,resU$wb.ratio,resU$sindex)
    names(res_vecI) = metrics
    names(res_vecU) = metrics
    df=rbind(df,res_vecI,res_vecU)
    rownames(df)[c(nrow(df)-1,nrow(df))]=c("Intrinsic","Unsupervised")
  }
  return(df)
}
#######################
## BRCA_pure results ##
#######################

# Table 3: Main BRCA_pure results
dir_name="BRCA_pure_med500_MAD50"

tab3=tabulateResults(dataset="BRCA_pure",med_filt=500,MAD_filt=50)
save(tab3,file="Results/Tab3.out")

####################################################################################################################

# Figure 3: BRCA_pure PAM50 heatmap

#read in results and data
load(sprintf("%s/FSC_covarsF_summary.out",dir_name)); resF=res
load(sprintf("%s/FSC_covarsT_summary.out",dir_name)); resT=res
load(sprintf("%s/FSC2_covarsT_summary.out",dir_name)); resT2=res
load("BRCA_pure_env.RData"); load("./TCGA_BRCA/BRCA_raw_normalizations.RData") # raw_norm_y plotted later to show all PAM50 genes

#filter raw_norm_y according to same criteria
idx=rowMeds>=500 & mads>=quantile(mads,50/100)                                 # gene pre-filtering criteria
colnames(raw_norm_y)=raw_anno$barcode; rownames(raw_norm_y)=rownames(raw_cts) # make sure column/row names are set
idy = colnames(raw_norm_y) %in% colnames(cts); raw_norm_y=raw_norm_y[,idy]     # match according purity filtered samples
colnames(norm_y)=colnames(cts)

#row annotations showing exclusions due to pre-filtering/FSCseq FS process
genelist = rep(NA,nrow(cts))
genelist_all = rep(NA,nrow(raw_cts))
for(i in 1:nrow(cts)){ genelist[i] = unlist(strsplit(rownames(cts)[i],"|",fixed=TRUE))[1] }
for(i in 1:nrow(raw_cts)){ genelist_all[i] = unlist(strsplit(rownames(raw_cts)[i],"|",fixed=TRUE))[1] }
pam50 = c("UBE2T","BIRC5","NUF2","CDC6","CCNB1","TYMS","MYBL2","CEP55","MELK","NDC80",
          "RRM2","UBE2C","CENPF","PTTG1","EXO1","ORC6L","ANLN","CCNE1","CDC20","MKI67",
          "KIF2C","ACTR3B","MYC","EGFR","KRT5","PHGDH","CDH3","MIA","KRT17","FOXC1",
          "SFRP1","KRT14","ESR1","SLC39A6","BAG1","MAPT","PGR","CXXC5","MLPH","BCL2",
          "MDM2","NAT1","FOXA1","BLVRA","MMP11","GPR160","FGFR4","GRB7","TMEM45B","ERBB2")
genelist2=genelist[idx]

PF_incl = as.factor(pam50 %in% genelist[idx]) # 41/50 passed pre-filtering
resF_incl = as.factor(pam50 %in% genelist2[resF$discriminatory])
resT_incl = as.factor(pam50 %in% genelist2[resT$discriminatory])
resT2_incl = as.factor(pam50 %in% genelist2[resT2$discriminatory])

annotation_row = data.frame(PF_incl,
                            resF_incl, resT_incl, resT2_incl)      # passed low count filtering (full), passed MAD pre-filtering, selected disc (disc)
colnames(annotation_row)=c("Passed PF","Disc (FSC)", "Disc (FSCadj)", "Disc (FSC10adj)")#expression(paste("Disc (",FSCadj^"10",")")))
rownames(annotation_row)=pam50
mycolors_PF=c("#FFFFFF","#000000")
mycolors_resF=c("#FFFFFF","#000000")
mycolors_resT=c("#FFFFFF","#000000")
mycolors_resT2=c("#FFFFFF","#000000")
names(mycolors_PF)=levels(PF_incl)
names(mycolors_resF)=levels(resF_incl)
names(mycolors_resT)=levels(resT_incl)
names(mycolors_resT2)=levels(resT2_incl)
mycolors2 = list("Passed PF"=mycolors_PF,
                 "Disc (FSC)"=mycolors_resF,
                 "Disc (FSCadj)"=mycolors_resT,
                 "Disc (FSC10adj)"=mycolors_resT2
                 #expression(paste("Disc (",FSCadj^"10",")"))=mycolors_resT2
                 )

# annotation col
# competing performers: KM, lMC
load("~/Research/P1/Real_Data/BRCA_pure3_2_med500_MAD50/KM_summary.out")
load("~/Research/P1/Real_Data/BRCA_pure3_2_med500_MAD50/lMC_summary.out")

annotation_col = data.frame(factor(lMC_summary$cls),
                            factor(KM_summary$cls),
                            factor(resF$clusters),
                            factor(resT$clusters),
                            factor(resT2$clusters),
                            factor(cls))
colnames(annotation_col)=c("lMC","KM","FSC","FSCadj","FSC10adj","Annotated")
rownames(annotation_col)=colnames(cts)

# order colors manually for highest consistency in coloring
newCols <- colorRampPalette(grDevices::rainbow(5))
mycolors_anno=newCols(3)
mycolors_FSC=newCols(5)[c(2,1,5,4,3)]
mycolors_FSCadj=newCols(3)[c(3,2,1)]
mycolors_FSC10adj=newCols(3)[c(3,2,1)]
mycolors_lMC=newCols(4)[c(2,1,4,3)]
mycolors_KM=newCols(2)[c(2,1)]
names(mycolors_anno)=unique(cls)[order(unique(cls))]
names(mycolors_FSC)=unique(resF$clusters)[order(unique(resF$clusters))]
names(mycolors_FSCadj)=unique(resT$clusters)[order(unique(resT$clusters))]
names(mycolors_FSC10adj)=unique(resT2$clusters)[order(unique(resT2$clusters))]
names(mycolors_lMC)=unique(lMC_summary$cls)[order(unique(lMC_summary$cls))]
names(mycolors_KM)=unique(KM_summary$cls)[order(unique(KM_summary$cls))]

mycolors1 = list(lMC=mycolors_lMC,KM=mycolors_KM,FSC=mycolors_FSC,FSCadj=mycolors_FSCadj,FSC10adj=mycolors_FSC10adj,Annotated=mycolors_anno)

# colors/breaks #
showpanel <- function(Colors){image(matrix(1:length(Colors), ncol=1), col=Colors, xaxt="n", yaxt="n" )}
oldpar <- par(mfrow=c(1,1))
my_cols<-colorRampPalette( c("yellow", "black", "cyan"), space="rgb")(100)
myBreaks <- seq(-2, 2, length.out=101)
dev.off()

# plot pheatmap PAM50 #
library(pheatmap)
rownames(raw_norm_y) = genelist_all
# png("Results/Fig3.png",height=850,width=800,res=300,type="cairo")
pheatmap(log(raw_norm_y[genelist_all %in% pam50,order(cls,resT$clusters,resF$clusters,lMC_summary$cls)]+0.1),scale="row",cluster_cols=F,
         annotation_col=annotation_col,annotation_row=annotation_row,annotation_colors=c(mycolors1,mycolors2), color=my_cols,show_colnames=F,
         breaks=myBreaks,main="TCGA BRCA High Purity Samples, PAM50 Genes with Annotated Exclusions from Pre-filtering and Analyses",
         height=12,width=11,filename="Results/Fig3.png")
# dev.off()

# for prelim presentation
# draw_colnames_45 <- function (coln, gaps, ...) {
#   coord <- pheatmap:::find_coordinates(length(coln), gaps)
#   x     <- coord$coord - 0.5 * coord$size
#   res   <- grid::textGrob(
#     coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
#     vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
#   )
#   return(res)
# }
# assignInNamespace(
#   x = "draw_colnames",
#   value = "draw_colnames_45",
#   ns = asNamespace("pheatmap")
# )
# png("C:/Users/limdd/Dropbox/Dissertation/Presentation/P1_PAM50.png",height=700,width=900)
# pheatmap(t(log(raw_norm_y[genelist_all %in% pam50,order(cls,resT$clusters,resF$clusters)]+0.1)),scale="column",cluster_rows=F,
#          annotation_row=annotation_col,annotation_col=annotation_row,annotation_colors=c(mycolors2,mycolors1), color=my_cols,show_rownames=F,
#          breaks=myBreaks,main="TCGA BRCA High Purity Samples, PAM50 Genes with Annotated Exclusions from Pre-filtering and Analyses",)
# dev.off()

####################################################################################################################

# Supp Figure 1: BRCA_pure all cluster-disc genes heatmap
load("./BRCA_pure_med500_MAD50/HC_summary.out")
load("./BRCA_pure_med500_MAD50/vMC_summary.out")
load("./BRCA_pure_med500_MAD50/rMC_summary.out")
load("./BRCA_pure_med500_MAD50/NBMB_summary.out")
load("./BRCA_pure_med500_MAD50/NMF_summary.out")
load("./BRCA_pure_med500_MAD50/iCl_summary.out")

annotation_col = data.frame(factor(rMC_summary$cls),
                            factor(vMC_summary$cls),
                            factor(lMC_summary$cls),
                            factor(NBMB_summary$cls),
                            factor(NMF_summary$cls),
                            factor(KM_summary$cls),
                            factor(HC_summary$cls),
                            factor(iCl_summary$cls),
                            factor(resF$clusters),
                            factor(resT$clusters),
                            factor(resT2$clusters),
                            factor(cls))
colnames(annotation_col)=c("rMC","vMC","lMC","NBMB","NMF","KM","HC","iCl","FSC","FSCadj","FSC10adj","Annotated")
rownames(annotation_col)=colnames(cts)

newCols <- colorRampPalette(grDevices::rainbow(8))
mycolors_anno=newCols(3)
mycolors_FSC=newCols(5)[c(2,1,5,4,3)]
mycolors_FSCadj=newCols(3)[c(2,3,1)]
mycolors_FSC10adj=newCols(3)[c(2,3,1)]
mycolors_lMC=newCols(4)[c(4,1,2,3)]
mycolors_rMC=newCols(4)[c(4,1,2,3)]
mycolors_vMC=newCols(4)[c(4,1,2,3)]
mycolors_KM=newCols(3)[c(2,1)]
mycolors_HC=newCols(3)[c(2,1)]
mycolors_iCl=newCols(8)[c(1:8)]
mycolors_NBMB=newCols(4)[c(1,2,4,3)]
mycolors_NMF=newCols(3)[c(2,1)]

names(mycolors_anno)=unique(cls)[order(unique(cls))]
names(mycolors_FSC)=unique(resF$clusters)[order(unique(resF$clusters))]
names(mycolors_FSCadj)=unique(resT$clusters)[order(unique(resT$clusters))]
names(mycolors_FSC10adj)=unique(resT2$clusters)[order(unique(resT2$clusters))]
names(mycolors_lMC)=unique(lMC_summary$cls)[order(unique(lMC_summary$cls))]
names(mycolors_vMC)=unique(vMC_summary$cls)[order(unique(vMC_summary$cls))]
names(mycolors_rMC)=unique(rMC_summary$cls)[order(unique(rMC_summary$cls))]
names(mycolors_KM)=unique(KM_summary$cls)[order(unique(KM_summary$cls))]
names(mycolors_HC)=unique(HC_summary$cls)[order(unique(HC_summary$cls))]
names(mycolors_iCl)=unique(iCl_summary$cls)[order(unique(iCl_summary$cls))]
names(mycolors_NBMB)=unique(NBMB_summary$cls)[order(unique(NBMB_summary$cls))]
names(mycolors_NMF)=unique(NMF_summary$cls)[order(unique(NMF_summary$cls))]


mycolors1 = list(rMC=mycolors_rMC,vMC=mycolors_vMC,lMC=mycolors_lMC,
                 NBMB=mycolors_NBMB,NMF=mycolors_NMF,KM=mycolors_KM,HC=mycolors_HC,
                 iCl=mycolors_iCl,FSC=mycolors_FSC,FSCadj=mycolors_FSCadj,FSC10adj=mycolors_FSC10adj,Annotated=mycolors_anno)

# plot pheatmap all disc genes from FSCadj #
# png("Results/SFig1.png",height=1000,width=800,res=300,type="cairo")
pheatmap(log(norm_y[idx,order(cls,resT$clusters,resF$clusters,lMC_summary$cls,vMC_summary$cls,rMC_summary$cls)][resT$discriminatory,]+0.1),scale="row",cluster_cols=F,
         annotation_col=annotation_col,annotation_colors=mycolors1, color=my_cols,show_colnames=F,show_rownames=F,
         breaks=myBreaks,main="TCGA BRCA High Purity Samples, Cluster-discriminatory Genes",height=14,width=11,filename="Results/SFig1.png")
# dev.off()

####################################################################################################################

# Supp Figure 2: TCGA Gene Ontology Analysis
library(TCGAbiolinks)
# Enrichment Analysis EA
# Gene Ontology (GO) and Pathway enrichment by DEGs list
Genelist <- genelist2[resF$discriminatory]
ansEA <- TCGAanalyze_EAcomplete(TFname="BRCA Disc genes found by FSC",Genelist)

TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = Genelist,
                        nBar = 10,
                        filename = "./Results/SFig2_FSC.pdf"
)
Genelist <- genelist2[resT$discriminatory]
ansEA <- TCGAanalyze_EAcomplete(TFname="BRCA Disc genes found by FSCadj",Genelist)

TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = Genelist,
                        nBar = 10,
                        filename = "./Results/SFig2_FSCadj.pdf"
)
####################################################################################################################

# Supp Figure 3: basalUP and basalDOWN, gene clustering --> GSEA Analysis
load("BRCA_pure_env.RData")
# cluster genes --> find genes that are upregulated/downregulated for basal subtype
library(MBCluster.Seq)
set.seed(999)
idx=rowMeds>=500 & mads>=quantile(mads,0.5)
mydata=RNASeq.Data(Count=norm_y[idx,][resT$discriminatory,],Normalize=NULL,Treatment=cls,GeneID=genelist[idx][resT$discriminatory])
c0=KmeansPlus.RNASeq(data=mydata,nK=10)$centers
fit=Cluster.RNASeq(data=mydata,centers=c0,model="nbinom",method="EM")

# manually determine which gene clusters are up/down regulated in basal subtype
pheatmap(log(norm_y[idx,][resT$discriminatory,][fit$cluster==6,order(cls)]+0.1),scale="row",cluster_cols=F,
         annotation_col=annotation_col,annotation_colors=mycolors1, color=my_cols,show_colnames=F,
         breaks=myBreaks,main="TCGA BRCA High Purity Samples, PAM50 Genes with Annotated Exclusions from Pre-filtering and Analyses")
basalUP_cls=c(5,6,10)
basalDOWN_cls=c(2,3,4,7,8)
lumADOWN_cls=c(1,9)

disc_genelistT=genelist2[resT$discriminatory]
disc_genelistF=genelist2[resF$discriminatory]
basalUP=disc_genelistT[fit$cluster%in% basalUP_cls]
basalDOWN=disc_genelistT[fit$cluster%in% basalDOWN_cls]


# load known gene sets to compare to
library(msigdbr)
m_df=msigdbr::msigdbr(species="Homo sapiens")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# Calculate p-values of overlap with each of the gene sets
library(GeneOverlap)
pvals_up=rep(NA,length(m_list))
pvals_down=rep(NA,length(m_list))
for(i in 1:length(m_list)){
  fitUP=newGeneOverlap(m_list[[i]],basalUP)
  fitUP=testGeneOverlap(fitUP)
  pvals_up[i]=fitUP@pval
  fitDOWN=newGeneOverlap(m_list[[i]],basalDOWN)
  fitDOWN=testGeneOverlap(fitDOWN)
  pvals_down[i]=fitDOWN@pval
}
pvals_up=p.adjust(pvals_up, method = "bonferroni", n = length(pvals_up))
pvals_down=p.adjust(pvals_down, method = "bonferroni", n = length(pvals_down))

head(names(m_list)[order(pvals_up)],n=10)
head(names(m_list)[order(pvals_down)],n=10)

# Plot top overlapping gene sets
library(ggplot2)
library(gridExtra)
bUP_df=data.frame(cbind(names(m_list)[order(pvals_up)],
                        as.numeric(-log(pvals_up[order(pvals_up)]))))
names(bUP_df)=c("Pathways","pval")
bUP_df2=bUP_df[1:10,]
bUP_df2$Pathways=factor(bUP_df2$Pathways,levels=bUP_df2$Pathways[10:1])
bUP_df2$pval=as.numeric(as.character(bUP_df2$pval))

p1=ggplot(bUP_df2, aes(x=Pathways, y=pval)) +
  geom_bar(stat = "identity") + ylab("-log(p-value)") + xlab("Pathways") +
  coord_flip()+ggtitle("basalUP")
# p1

bDOWN_df=data.frame(cbind(names(m_list)[order(pvals_down)],
                          as.numeric(-log(pvals_down[order(pvals_down)]))))
names(bDOWN_df)=c("Pathways","pval")
bDOWN_df2=bDOWN_df[1:10,]
bDOWN_df2$Pathways=factor(bDOWN_df2$Pathways,levels=bDOWN_df2$Pathways[10:1])
bDOWN_df2$pval=as.numeric(as.character(bDOWN_df2$pval))

p2=ggplot(bDOWN_df2, aes(x=Pathways, y=pval)) +
  geom_bar(stat = "identity") + ylab("-log(p-value)") + xlab("Pathways") +
  coord_flip()+ggtitle("basalDOWN")
# p2

# png("./Results/SFig3.png",width=900,height=900,res=300,type="cairo")
# grid.arrange(p1,p2)
# dev.off()

p = arrangeGrob(p1,p2,nrow=2)
ggsave(file="./Results/SFig3.png",width=10,height=10, p)

####################################################################################################################

#######################
## BRCA_full results ##
#######################

dir_name="BRCA_full_med500_MAD50"

# Supp Table 8: BRCA_full results table
stab8=tabulateResults(dataset="BRCA_full",med_filt=500,MAD_filt=50)
save(stab8,file="Results/STab8.out")

####################################################################################################################

# Supp Figure 4: BRCA_full all cluster-disc genes heatmap

load("./BRCA_full_med500_MAD50/FSC_covarsF_summary.out"); resF=res
load("./BRCA_full_med500_MAD50/FSC2_covarsF_summary.out"); resF2=res
load("./BRCA_full_med500_MAD50/FSC_covarsT_summary.out"); resT=res
load("./BRCA_full_med500_MAD50/FSC2_covarsT_summary.out"); resT2=res
load("BRCA_full_env.RData")
idx=rowMeds>=500 & mads>=quantile(mads,0.5)
# annotation col
# competing best performer: KM
load("./BRCA_full_med500_MAD50/KM_summary.out")

annotation_col = data.frame(factor(anno$cls_sigCL_intrinsic),
                            factor(anno$cls_sigCL_unsupervised),
                            factor(KM_summary$cls),
                            factor(resF$clusters),
                            factor(resF2$clusters),
                            factor(resT$clusters),
                            factor(resT2$clusters),
                            factor(cls))
colnames(annotation_col)=c("SigI","SigU",
                           "KM","FSC","FSC10","FSCadj","FSC10adj","Anno")
rownames(annotation_col)=colnames(cts)

newCols <- colorRampPalette(grDevices::rainbow(14))
mycolors_anno=newCols(5)
mycolors_FSC=newCols(7)[c(4,2,3,1,5,6,7)]
mycolors_FSC10=newCols(5)[c(4,5,1,3,2)]
mycolors_FSCadj=newCols(8)[c(3,2,5,7,4,6,1,8)]
mycolors_FSC10adj=newCols(5)[c(3,4,5,1,2)]
mycolors_I=newCols(14)[c(5,2,3,4,1,6:14)]
mycolors_U=newCols(13)
mycolors_KM=newCols(14)[c(12,1)]
names(mycolors_anno)=unique(cls)[order(unique(cls))]
names(mycolors_FSC)=unique(resF$clusters)[order(unique(resF$clusters))]
names(mycolors_FSC10)=unique(resF2$clusters)[order(unique(resF2$clusters))]
names(mycolors_FSCadj)=unique(resT$clusters)[order(unique(resT$clusters))]
names(mycolors_FSC10adj)=unique(resT2$clusters)[order(unique(resT2$clusters))]
names(mycolors_I)=unique(anno$cls_sigCL_intrinsic)[order(unique(anno$cls_sigCL_intrinsic))]
names(mycolors_U)=unique(anno$cls_sigCL_unsupervised)[order(unique(anno$cls_sigCL_unsupervised))]
names(mycolors_KM)=unique(KM_summary$cls)[order(unique(KM_summary$cls))]

mycolors1 = list(SigI=mycolors_I, SigU=mycolors_U,
                 KM=mycolors_KM,FSC=mycolors_FSC,FSC10=mycolors_FSC10,FSCadj=mycolors_FSCadj,FSC10adj=mycolors_FSC10adj,Anno=mycolors_anno)

# colors/breaks #
showpanel <- function(Colors){image(matrix(1:length(Colors), ncol=1), col=Colors, xaxt="n", yaxt="n" )}
oldpar <- par(mfrow=c(1,1))
my_cols<-colorRampPalette( c("yellow", "black", "cyan"), space="rgb")(100)
myBreaks <- seq(-2, 2, length.out=101)
dev.off()

# png("./Results/SFig4.png",height=980,width=930,res=300,type="cairo")
pheatmap(log(norm_y[idx,order(cls,resT2$clusters,resT$clusters)][resT$discriminatory,]+0.1),scale="row",cluster_cols=F,
         annotation_col=annotation_col,annotation_colors=mycolors1, color=my_cols,show_colnames=F,show_rownames=F,
         breaks=myBreaks,main="TCGA BRCA All Samples, Disc Genes from FSCseq, Ordered by Subtype",height=14,width=12,filename="./Results/SFig4.png")
# dev.off()


####### Correlation Plots #######
corPlots=function(dataset=c("BRCA_pure","BRCA_plate"),covars){
  library(corrplot)
  load(sprintf("%s_env.RData",dataset))
  filt_idx = rowMeds>=500 & mads>=quantile(mads,0.5)
  norm_y=norm_y[filt_idx,]
  load(sprintf("%s_med500_MAD50/FSC_covars%s_summary.out", dataset,covars))
  discriminatory=res$discriminatory
  disc_norm_y = norm_y[discriminatory,]
  cors = cor(t(disc_norm_y))
  covar_text1 = if(dataset=="BRCA_pure"){"BRCA (high purity samples)"} else{"BRCA (all samples)"}
  covar_text2 = if(covars=="T"){"Adjusted for Plate"} else{"Not adjusted for Plate"}

  png(sprintf("%s_med500_MAD50/CEM_covars%s_corPlot.png",dataset,covars),res=300,width=1200,height=1300,type="cairo",pointsize=6); par(xpd=TRUE)
  corrplot(cors,method="color",tl.pos = "td",tl.cex = 0.2,mar=c(0,0,3,0),tl.col = 'black',type="upper",order="hclust",diag=F,
           main=sprintf("Correlations of %s cluster-discriminatory genes \n%s",covar_text1,covar_text2))
  dev.off()
}
corPlots("BRCA_pure","F")
corPlots("BRCA_pure","T")
corPlots("BRCA_full","F")
corPlots("BRCA_full","T")

####### FSC2 and FSCadj2 #######
# finding smallest model in bottom 10% of BIC
opt_model = function(dataset,covars){
  K_search=2:8; lambda_search=seq(0.25,5,0.25); alpha_search=c(0.01,seq(0.05,0.5,0.05))
  npoints = length(K_search)*length(lambda_search)*length(alpha_search)
  BICs=matrix(NA,nrow=npoints,ncol=6)
  index=1
  colnames(BICs) = c("K","l","a","BIC","ndisc","nparams")
  for(c in 1:length(K_search)){for(l in 1:length(lambda_search)){for(a in 1:length(alpha_search)){
    BICs[index,1:3] = c(K_search[c], lambda_search[l],alpha_search[a])
    file.name=sprintf("%s_med500_MAD50/joint%d_%f_%f_gene_CEM_covars%s.out", dataset,K_search[c],lambda_search[l],alpha_search[a],covars)
    if(file.exists(file.name)){load(file.name); BICs[index,4:6] = c(res$BIC, sum(res$discriminatory), res$num_est_params)}else{print(file.name)}
    index=index+1
  }}}

  #BICs[order(BICs[,4])[1:100,],]

  return(
    BICs[order(BICs[,4]),] [which.min(BICs[order(BICs[,4])[1:floor(npoints/10)],][,6]),]    # lowest nparams from bottom 10% of BIC
  )
}
p1=opt_model("BRCA_pure","F")
p2=opt_model("BRCA_pure","T")
# covarF : 2.00       0.25       0.50 8704580.49    1169.00
# covarT : 3          1.00       0.30 9789277       1324
f1=opt_model("BRCA_full","F")
f2=opt_model("BRCA_full","T")
# covarF: 5.0         1.5         0.1  42533187.3      2800.0
# covarT: 5.00        3.25        0.05 43079204.94     2834.00

check_opt_model = function(dataset,covars,K,lambda,alpha){
  load(sprintf("%s_med500_MAD50/joint%d_%f_%f_gene_CEM_covars%s.out",dataset,K,lambda,alpha,covars))
  load(sprintf("%s_env.RData",dataset))
  library(mclust)
  return(list(ARI=adjustedRandIndex(res$clusters,cls),
              ndisc=sum(res$discriminatory)))
}
check_opt_model("BRCA_pure","F",p1[1],p1[2],p1[3])
check_opt_model("BRCA_pure","T",p2[1],p2[2],p2[3])
# covarF : 2.00       0.25       0.50 8704580.49    1169.00   9246.00    ----> ARI = 0.486
# covarT : 3          1.00       0.30 9789277       1324      106377.0   ----> ARI = 0.610
check_opt_model("BRCA_full","F",f1[1],f1[2],f1[3])
check_opt_model("BRCA_full","T",f2[1],f2[2],f2[3])
# covarF: 5.0        1.5        0.1 42533187.3     2800.0    12867.0  ----> ARI = 0.254
# covarT: 5.0        1.5        0.1 43072187.8     2900.0    59580.0  ----> ARI = 0.288
