
# Real Data

FSCseq analysis was performed on the TCGA Breast Cancer dataset

## Acquiring Data

Run `./TCGA_BRCA/DataAcquisition.R` to acquire the dataset and
annotations. The dataset used in our analyses was downloaded on March
25th, 2019.

## Pre-processing Data

Run `./TCGA_BRCA/DataPreProcessing.R` to pre-process the data in the
same way as in the manuscript. Cluster labels `SigI` and `SigU` that
were compared in the FSCseq paper were extracted from annotations from
the TCGA BRCA [paper](https://doi.org/10.1038/nature11412). Batch and
plate information was downloaded manually from
[here](https://bioinformatics.mdanderson.org/BatchEffectsViewer/), with
Index File `GDC 2019-07-30-1200 (current)`, and Workflow
`RNAseq-counts`. This information was saved to
`./RealData/TCGA_BRCA/BatchData.tsv`. Purity information of samples were
extracted from this [paper](https://doi.org/10.1038/ncomms9971) by Aran
et al (2014). In FSCseq, we analyze `BRCA_full` and `BRCA_pure`, and
necessary objects after pre-filtering, including normalized counts, are
saved in `BRCA_full_env.RData` and `BRCA_pure_env.RData`, respectively,
for ease of use in submitting slurm batch jobs. Data without applying
pre-filtering steps is also saved in
`./RealData/TCGA_BRCA/BRCA_raw_normalizations.RData` for use in creating
some post-analysis figures.

## Creating R Scripts and Submitting Batch Jobs via Slurm

R scripts for all jobs pertaining to real data are created and submitted
as batch jobs to slurm using the `runRealAnalyses_submitScripts_*.R`
files, which sources functions from `runRealAnalyses_runScripts.R`. For
computational efficiency, the FSCseq workflow was broken down into two
sequential parts: 1. `*FSCinit.R` submits jobs for the initialization
searches with small penalty, and 2. `*FSCtune.R` searches across values
of tuning parameters ![K](https://latex.codecogs.com/png.latex?K "K"),
![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\\lambda"),
and ![\\alpha](https://latex.codecogs.com/png.latex?%5Calpha "\\alpha").
The `*others.R` submits jobs for competing methods.

## Collecting Results

Results are saved in subdirectories corresponding to the dataset name,
and gene pre-filtering thresholds for median count and MAD quantile,
which are fixed in FSCseq for all real data analyses at
![500](https://latex.codecogs.com/png.latex?500 "500") and
![50](https://latex.codecogs.com/png.latex?50 "50"), respectively. Run
`collectFSCResults.R` to compare the BICs from all FSCseq tuning
results, and select the optimal values for tuning parameters. For easy
access in summary of results, copy the corresponding result into a new
file `FSC_covars[*]_summary.out`, with `[*]` denoting whether the plate
effect was adjusted.

## Leave-one-out Cross Validation

After `FSC_res_covars[*].out` are saved, run leave-one-out
cross-validation using `Leave1Out_CV.R`.

## Replicating Tables and Figures

Run `SummarizeResults.R` to replicate tables and figures pertaining to
`BRCA_full` and `BRCA_pure` in the FSCseq paper. Output is directed into
the `/Results` subdirectory.
