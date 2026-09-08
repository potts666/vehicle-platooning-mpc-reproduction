# Vehicle Platooning MPC Reproduction

Small, self-contained MATLAB implementations for three vehicle-platooning MPC studies. This repository contains source code and run instructions only; generated MATLAB data, figures, papers, reading notes, reports, caches, local environments, and credentials are deliberately excluded.

## Contents

| Directory | Entry point | Scope |
| --- | --- | --- |
| `zheng-2017-recursive-feasibility/` | `step1_nominal_mpc.m` through `step9_constant_spacing_rmpc_vs_dmpc.m` | A staged study of nominal MPC, terminal sets, robust MPC, and communication-based DMPC. |
| `zheng-2017-heterogeneous-topologies/` | `zheng_pf_nonlinear_dmpc.m` | A nonlinear heterogeneous-platoon simulation with selectable unidirectional communication topology. |
| `sun-2026-asynchronous-communication/` | `sun_2026_abd_dmpc_validation.m` | A carefully scoped implementation and numerical audit for asynchronous communication with coupled constraints. |

## Prerequisites

See [environment.md](environment.md). MATLAB scripts are run from their own project directory so that generated `.mat` files remain local and are ignored by Git.

## Reproducibility boundaries

These scripts are research reproductions and engineering studies, not a substitute for the cited papers or their formal proofs. In particular, the Sun project distinguishes a feasible engineering experiment from literal-paper numerical auditing and from an exploratory extension. Its README states those boundaries explicitly.

## Repository hygiene

Never commit generated results, local configuration, credentials, or papers/reports to this repository. The `.gitignore` covers the common artifacts; review `git status` before every commit.
