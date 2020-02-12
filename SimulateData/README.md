
# Simulate Data

Run `SimulateRNASeqData.R` to simulate datasets spanning all of the
conditions in the FSCseq paper. Simulated data will be saved in a
subdirectory in “./Simulations/”, which will be named according to the
simulated parameters. For example, a dataset with parameters
{`sigma_g=0.1`, `sigma_b=0`, `B=1`, `LFCb=0`, `K=2`, `n=50`, `LFCg=1`,
`pDEg=0.025`, `beta0=8`, and `phi0=0.15`} with `sim_index=1` will be
saved
as:

`./Simulations/0.100000_0.000000/B1_LFCb0/2_50_1.000000_0.025000_8.000000_0.150000_sim1_data.RData`

Subdirectories that don’t exist will be automatically created.
