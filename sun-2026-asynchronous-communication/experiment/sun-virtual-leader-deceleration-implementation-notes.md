# Virtual-leader deceleration: implementation notes

- Scenario: virtual leader jumps from 10 m/s to 10.5 m/s at t = 0,
  then decelerates uniformly to 9.5 m/s over 2 s.
- Initial run: stopped before optimization because `eq_torque` received
  a vector-valued speed reference.
- Cause: the expression used `v^2`, which only supports scalar `v`.
- Resolution: changed it to element-wise `v.^2`.
- Interpretation: this was an implementation correction, not a
  numerical-feasibility result.