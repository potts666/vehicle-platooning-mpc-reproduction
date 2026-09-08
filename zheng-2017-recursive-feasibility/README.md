# Zheng 2017: recursive-feasibility study

This directory is a staged MATLAB study of nominal MPC, terminal sets, robust MPC, and communication-based DMPC. It is intended to be run from this directory.

## Requirements

MATLAB with Control System Toolbox and Optimization Toolbox, plus YALMIP and MPT3. Initialize MPT3 (`mpt_init`) and ensure YALMIP is on the MATLAB path before running the scripts.

## Run order

Run the scripts in numerical order:

1. `step1_nominal_mpc.m`
2. `step2_terminal_set.m`
3. `step3_mpc_with_terminal_set.m`
4. `step4_robust_invariant_set.m`
5. `step5_pontryagin_difference.m`
6. `step6_rmpc_with_worst_case_disturbance.m`
7. `step7_dmpc_with_communicated_input.m`
8. `step8_compare_feasible_regions.m`
9. `step9_constant_spacing_rmpc_vs_dmpc.m`

Steps 3–9 load `.mat` outputs created by earlier steps. Those files are deliberately **not** in the repository: run the prerequisite steps locally to regenerate them. In particular, Steps 3, 6, 7, and 8 require `terminal_set_step2.mat`; Step 5 requires `robust_sets_step4.mat`; Steps 6 and 8 require `rmpc_constraint_step5.mat`; Step 7 can use Step 6's simulation output for comparison.

The scripts produce figures and additional `.mat` artifacts locally. They are ignored by Git.
