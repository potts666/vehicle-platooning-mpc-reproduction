# Environment

## Tested baseline

- MATLAB R2022a on Windows.
- MATLAB Optimization Toolbox (`quadprog` and `fmincon`).
- MATLAB Control System Toolbox (`dlqr` and `ctrb`) for the Zheng 2017 recursive-feasibility steps.
- YALMIP and MPT3 for the Zheng 2017 recursive-feasibility steps.

## Project-specific requirements

| Project | MATLAB dependencies |
| --- | --- |
| `zheng-2017-recursive-feasibility` | Control System Toolbox, Optimization Toolbox, YALMIP, MPT3 |
| `zheng-2017-heterogeneous-topologies` | Optimization Toolbox (`fmincon`) |
| `sun-2026-asynchronous-communication` | Optimization Toolbox (`fmincon`, `lsqnonlin`) |

Before running the staged Zheng project, add YALMIP and MPT3 to the MATLAB path and run `yalmip('clear')`/`mpt_init` as appropriate for the installed versions. The code itself does not install or configure dependencies.

Generated `.mat`, `.fig`, logs, and Sun result files are local execution artifacts and are intentionally excluded from version control.
