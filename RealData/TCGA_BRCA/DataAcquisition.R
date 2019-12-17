# library(devtools)
# devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
TCGAbiolinks:::getGDCprojects()$project_id
query1=GDCquery(project="TCGA-BRCA",
                data.category = "Gene expression",
                data.type = "Gene expression quantification",
                platform = "Illumina HiSeq",
                file.type  = "results",
                experimental.strategy = "RNA-Seq",
                legacy = TRUE)
GDCdownload(query1)
GDCprepare(query = query1,
           save = TRUE,
           save.filename = "TCGA_BRCA_exp.rda")

query2=GDCquery(project="TCGA-BRCA",
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification",
                workflow.type = "HTSeq - FPKM")
BRCA_FPKM = getResults(query2)

GDCdownload(query=query2,directory="C:/Users/limdd/Documents/Research/")
GDCprepare(query=query2,
           save = TRUE, directory="C:/Users/limdd/Documents/Research/",
           save.filename="TCGA_BRCA_FPKM.rda",summarizedExperiment = TRUE)
