
# Simulations

## Simulating Data

See
[`SimulateRNASeqData`](https://github.com/DavidKLim/FSCseqPaper/tree/master/SimulateData).
Simulated datasets will be saved in a subdirectories of “./Simulations”

## Creating R Scripts and Submitting Batch Jobs via Slurm

Run `submitJobs_slurm.R` to both create scripts (using
`runSimAnalyses_createScripts.R`, saved to “scripts/” directory) and
submit R batch jobs into slurm (`runSimAnalyses_runScripts_slurm.R`).
Jobs whose results have already been attained will be skipped. Each
created script will use wrapper functions from `runSimAnalyses.R`, which
use the functions in `runSimAnalyses_functions.R` to perform respective
analyses.

## Collecting Results

Results from analyses will be saved in the same subdirectories that
contain the simulated data. Run `collectSimResults.R` to collect
results. This will save results into the “./Simulations/Results/”
directory, and we will use those files to produce the tables and figures
in the text.

## Replicating Tables and Figures

After all results are finished and saved, run `SummarizeResults.R` to
replicate tables and figures from simulations in the FSCseq paper. Ouput
is directed to the `/Results`
[subdirectory](https://github.com/DavidKLim/FSCseqPaper/tree/master/Simulations/Results)
