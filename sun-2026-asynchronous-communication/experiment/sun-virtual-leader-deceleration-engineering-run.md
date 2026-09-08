# Virtual-leader deceleration: engineering run

## Scenario

- Mode: `engineering`.
- Simulation duration: 10 s.
- Before `t = 0`, the virtual leader speed is 10 m/s.
- At `t = 0`, it jumps to 10.5 m/s.
- Over `0 <= t < 2 s`, it decelerates uniformly at -0.5 m/s^2.
- From `t = 2 s` onward, it remains at 9.5 m/s.

## Implementation note

The time-varying speed reference made `eq_torque` receive a vector. The
drag term was corrected from `v^2` to `v.^2` before this run so torque
references are evaluated element-wise for the full prediction horizon.

## Recorded outcome

The simulation completed to 10 s without an enforced physical-bound or gap
violation. The reported final tracking errors were
`[-1.7739e-07, 7.3464e-10, -2.3394e-07, 1.35e-13] m`.

| Metric | Recorded value |
| --- | ---: |
| Peak absolute tracking errors, vehicles 1--4 | 0.18185, 0.18903, 0.18725, 0.19010 m |
| Actual adjacent-gap range | 9.993938--10.181848 m |
| Allowed gap range | 5--15 m |
| Maximum accepted scaled constraint residual | 2.77e-08 |
| Script acceptance tolerance | 2e-06 |
| Feasible-fallback local updates | 143 |
| Vehicle 3 fallback interval/count | 2.3--9.9 s / 77 |
| Vehicle 1 fallback interval/count | 3.4--9.9 s / 66 |

## Interpretation and limits

At 143 local updates, `fmincon` returned `exitflag = -2`. The wrapper did
not apply that unsuccessful optimizer output. Instead, it applied the
independently checked feasible seed and issued a `SUN:FeasibleFallback`
warning. Therefore, this is an engineering simulation with verified-feasible
fallback controls, not a run in which every local problem returned an
accepted optimum.

The measured follower peak ordering was `0`; it does not demonstrate string
stability. This `engineering` result is neither a complete paper
reproduction nor a theoretical proof.

Generated `.mat`, `.fig`, and `.png` artifacts remain local under
`sun_abd_results/` and are intentionally not committed.
