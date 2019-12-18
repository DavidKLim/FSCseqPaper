library(SummarizedExperiment)
library(matrixStats)
# load downloaded data
load(file="./RealData/TCGA_BRCA/TCGA_BRCA_FPKM.rda")
FPKMdat=data; rm(data)
load(file="./RealData/TCGA_BRCA/TCGA_BRCA_exp.rda")

raw_cts <- round(assay(data),0)
anno <- colData(data)@listData
raw_cts <- raw_cts[!duplicated(raw_cts[,1:ncol(raw_cts)]),]    # remove duplicated genes

# Add annotated subtypes, fewer missing
library(TCGAbiolinks)
BRCA_tab = TCGAquery_subtype("BRCA")
anno$subtypes=rep(NA,ncol(raw_cts))
match_ids=match(anno$patient,BRCA_tab$patient)   # 1215: places
anno$subtypes=BRCA_tab$BRCA_Subtype_PAM50[match_ids]

# Added purity information from Aran et al (2015)
library(openxlsx)
tab_purity = read.xlsx("./RealData/TCGA_BRCA/TCGA_purity.xlsx",startRow=4)
tab_purity = tab_purity[tab_purity$Cancer.type=="BRCA",]
match_ids=match(anno$sample,tab_purity$Sample.ID)
anno$purity=as.numeric(tab_purity$ESTIMATE[match_ids])

# Read SigI and SigU from supplementary material of TCGA BRCA paper (2012)
library(readxl)
anno_paper=readxl::read_xls("./RealData/TCGA_BRCA/Supplementary Tables 1-4.xls",sheet=1,skip=1)
match_ids=match(anno$patient,anno_paper$`Complete TCGA ID`)
anno$cls_sigCL_unsupervised = anno_paper$`SigClust Unsupervised mRNA`[match_ids]
anno$cls_sigCL_intrinsic = anno_paper$`SigClust Intrinsic mRNA`[match_ids]
anno$cls_PAM50 = anno_paper$`PAM50 mRNA`[match_ids]

## FPKM pre-filtering by median(FPKM)>1 (match by gene name)
FPKM_cts = assay(FPKMdat)
rownames(FPKM_cts) = rowData(FPKMdat)$external_gene_name
genelist = rep(NA,nrow(raw_cts))
for(i in 1:nrow(raw_cts)){
  genelist[i] = unlist(strsplit(rownames(raw_cts)[i],"|",fixed=TRUE))[1]
}
#sum(rowData(FPKMdat)$external_gene_name %in% genelist)          # 17272 genes with FPKM data
#sum(rownames(colData(FPKMdat)) %in% anno$barcode)               # 1205 samples with FPKM data
match_idx = match(genelist,rownames(FPKM_cts))               # 2424 missing
match_idy = match(anno$barcode,colnames(FPKM_cts))           # 10 missing
FPKM_cts = FPKM_cts[match_idx, match_idy]
FPKM_count_filt_ids = rowMedians(FPKM_cts,na.rm=T)>1     # T/F/NA. Leave NA's, T --> stays, F --> filtered out
FPKM_idx = is.na(FPKM_count_filt_ids) | FPKM_count_filt_ids

# Add batch/plate information
BRCA_batch = read.table("./RealData/TCGA_BRCA/BatchData.tsv",header=TRUE,sep="\t")
match_id = match(anno$barcode,BRCA_batch$Sample)
anno$batch=BRCA_batch$BatchId[match_id]
anno$plate=BRCA_batch$PlateId[match_id]

## Raw normalizations (for visualization of PAM50 later) ##
# library(FSCseq)
# fit <- processData(raw_cts)
# raw_SF <- fit$size_factors
# raw_norm_y <- fit$norm_y
# raw_anno=anno
# save(list=c("raw_anno", "raw_cts", "raw_norm_y", "raw_SF"),file="./RealData/TCGA_BRCA/BRCA_raw_normalizations.RData")
load("./RealData/TCGA_BRCA/BRCA_raw_normalizations.RData")


###############
## BRCA_pure ##
###############
idy=!is.na(raw_anno$purity) & raw_anno$purity>0.9  # 134 samples
table(raw_anno$subtypes[idy]) # low incidence of Her2 (n=6) and Normal (n=3)

idy=idy & raw_anno$subtypes%in% c("Basal","LumA","LumB")
table(raw_anno$subtypes[idy]) # Filtered out low incidence subtypes

cts=raw_cts[,idy]; norm_y=raw_norm_y[,idy]; SF=raw_SF[idy]

cts=cts[FPKM_idx,]
norm_y=norm_y[FPKM_idx,]

anno=lapply(raw_anno,function(x) x[idy])
anno$plate = as.character(droplevels(anno$plate))
anno$plate[is.na(anno$plate)]=substr(anno$barcode[is.na(anno$plate)],22,25)

# merge plates with 1 sample
singleton_plates = names(table(anno$plate))[table(anno$plate)==1]
anno$plate[anno$plate %in% singleton_plates]="ones"

X_plate=matrix(0,nrow=ncol(cts),ncol=length(unique(anno$plate))-1)
for(i in 1:ncol(X_plate)){
  X_plate[,i] = (anno$plate==unique(anno$plate)[i+1])^2
}

X=X_plate

rowMeds=rowMedians(norm_y)
mads=rep(0,nrow(norm_y))
for(j in 1:nrow(norm_y)){
  mads[j] = mad(log(norm_y[j,]+0.1))
}

cls=anno$subtypes

save(list=c("cts","norm_y","SF","mads","rowMeds","anno","X","cls"),file="./RealData/BRCA_pure_env.RData")

###############
## BRCA_full ##
###############
# filter out just samples with batch/plate/subtype information that were included in the TCGA study
idy = !(is.na(raw_anno$batch) |
          is.na(raw_anno$plate) |
          is.na(raw_anno$cls_PAM50) | raw_anno$cls_PAM50 =="NA" |
          is.na(raw_anno$cls_sigCL_intrinsic) | raw_anno$cls_sigCL_intrinsic =="NA" |
          is.na(raw_anno$cls_sigCL_unsupervised) | raw_anno$cls_sigCL_unsupervised =="NA" |
          is.na(raw_anno$subtypes))

cts=raw_cts[,idy]; norm_y=raw_norm_y[,idy]; SF=raw_SF[idy]

anno=lapply(raw_anno,function(x) x[idy])
anno$plate = as.character(droplevels(anno$plate))

# agglomerate plates with 2 or fewer samples
doubleton_plates = names(table(anno$plate))[table(anno$plate)<=2]
anno$plate[anno$plate %in% doubleton_plates]="doubleton"

# cts=cts[,idy2]; norm_y=norm_y[,idy2]; SF=SF[idy2]
# anno3=lapply(anno2,function(x) x[idy2])
# anno3$plate = droplevels(anno3$plate); anno3$batch=droplevels(anno3$batch)


X_plate=matrix(0,nrow=ncol(cts),ncol=length(unique(anno$plate))-1)
for(i in 1:ncol(X_plate)){
  X_plate[,i] = (anno$plate==unique(anno$plate)[i+1])^2
}

cts=cts[FPKM_idx,]
norm_y=norm_y[FPKM_idx,]

rowMeds=rowMedians(norm_y)
mads=rep(0,nrow(norm_y))
for(j in 1:nrow(norm_y)){
  mads[j] = mad(log(norm_y[j,]+0.1))
}

cls=anno$subtypes
X=X_plate

save(list=c("cts","norm_y","SF","mads","rowMeds","anno","X","cls"),file="./RealData/BRCA_full_env.RData")

