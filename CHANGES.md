CHANGES
=======

Scope
-----
Delta from the original repo state to the current stable filter that passes
`tester.sh` consistently.

State and model
---------------
Original:
- State x = [p, v, a, dir] (12D).
- Direction was a measured state; acceleration and direction both updated as
  measurements.
- In BONUS, GPS noise was adaptively scaled (Huber/alpha logic).

Current:
- State x = [p, v, a] (9D).
- Direction is used only to project/align longitudinal motion; it is not a
  state or a measurement update.
- Motion is a constant-acceleration model driven by IMU updates.

Prediction (unchanged form, tightened to one loop)
--------------------------------------------------
x_k^- = F x_{k-1}
P_k^- = F P_{k-1} F^T + Q

Q uses the continuous white-noise acceleration integration:
q11 = q dt^5/20, q12 = q dt^4/8, q13 = q dt^3/6,
q22 = q dt^3/3, q23 = q dt^2/2, q33 = q dt

Measurement update (made explicit and stable)
---------------------------------------------
y = z - H x^-
S = H P^- H^T + R_eff
K = P^- H^T S^{-1}
x = x^- + K y
P = (I - K H) P^- (I - K H)^T + K R_eff K^T   (Joseph form)

GPS robustness (now spec-tied and minimal)
------------------------------------------
Original:
- BONUS used alpha-up/down and Huber gating.
- Base R_gps effectively huge due to n^7 scaling, which reduced GPS influence.

Current:
- scale = max(1, (||y||/t)^2, (NIS/dof) * (5/t))
- first GPS update enforces scale >= (NIS/dof)^2 (fast initial skepticism)
- R_eff = R_gps * scale
- Base R_gps also scaled by (5/t)^2 so tolerance directly controls trust.

Ordering and sensor usage
-------------------------
Original:
- Predict and update were split between `main.cpp` and `parse`, and the order
  varied across refactors.

Current:
- Predict then update within `parse`.
- Acceleration updates occur every IMU message.
- GPS update occurs when present.
- Direction only aligns velocity/acceleration longitudinally.

Diagnostics
-----------
- Added message-count logging to correlate failures with early GPS bursts.
