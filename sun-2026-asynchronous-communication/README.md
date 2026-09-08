# Sun 2026: asynchronous communication and coupled constraints

Run `sun_2026_abd_dmpc_validation.m` from this directory in MATLAB with Optimization Toolbox (`fmincon`, `lsqnonlin`). Generated diagnostics and results are written to `sun_abd_results/` and are intentionally ignored by Git.

## Modes and claims

The script has three explicitly different modes selected by `SUN_ABD_MODE` (or by editing `mode` near the top of the script):

- `engineering` (default) is a reproducible **engineering simulation** with staggered nonlinear MPC, bounds, collision constraints, and a hard terminal equality. It does **not** enforce the paper's string constraints and does **not** establish the paper's string-stability theorem.
- `paper_audit` applies a literal numerical audit of the paper constraints, including the first-follower condition. Under the current initialization and stated assumptions, it has **not** found a solution meeting the acceptance tolerance. It must not be reported as a numerically feasible replication.
- `anchored_string` is an **exploratory extension**: it introduces an explicit first-follower error envelope and begins the predecessor string chain at vehicle 2. It is not a complete reproduction of the paper and provides no theoretical proof.

Therefore, this project must not be described as a complete reproduction or as a theoretical verification of the paper.

## Run controls

- `SUN_ABD_TSIM`: positive simulation duration, such as `10`.
- `SUN_ABD_MODE`: `engineering`, `paper_audit`, or `anchored_string`.
- `SUN_ABD_NO_PLOTS`: any nonempty value suppresses plots.
- `SUN_ABD_E1MAX`: positive first-follower envelope for `anchored_string`.

Example (MATLAB):

```matlab
setenv('SUN_ABD_MODE', 'engineering');
setenv('SUN_ABD_TSIM', '10');
run('sun_2026_abd_dmpc_validation.m')
```
