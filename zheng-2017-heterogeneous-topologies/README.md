# Zheng 2017: heterogeneous vehicle platoons under unidirectional topologies

Run `zheng_pf_nonlinear_dmpc.m` from this directory in MATLAB with Optimization Toolbox installed.

The script simulates one leader and seven heterogeneous followers using nonlinear MPC. Its default topology is `TPLF`. Change the `topology` setting near the top of the script to select `PF`, `PLF`, `TPF`, or `TPLF`.

For a short check, set `TsimOverride` in the MATLAB workspace before running the script, for example `TsimOverride = 0.5;`. Generated figures and workspace artifacts are not versioned.
