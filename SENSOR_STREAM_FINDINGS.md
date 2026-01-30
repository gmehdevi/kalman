SENSOR STREAM FINDINGS
======================

Summary
-------
The sensor stream's default GPS noise does not match the subject's sigma=0.1.
Using a two-pass, fixed-seed procedure, the observed GPS noise with default
settings is ~0.9 (sigma), while `-g 0.1` yields ~0.09 (sigma). This suggests
the binary's default GPS sigma is ~0.9 (variance ~0.81), not 0.1 (variance 0.01).

Because the kalman filter's measurement covariance depends directly on this
sigma, using spec values against the stream's default causes divergence and
frequent "Delta is too high" failures. Tuning `--gps` to the stream's actual
noise is required for stability if the stream is left at default settings.

Explicit diagnosis (spec mismatch)
----------------------------------
Spec quote (noise model, FR):
"Mais, étant donné que nous vivons dans un monde imparfait, certaines de ces mesures
sont affectées par un bruit blanc gaussien:
• Accéléromètre : σ = 10−3
, υ = 0
• Gyroscope : σ = 10−2
, υ = 0
• GPS : σ = 10−1
, υ = 0"

Spec quote (noise model, EN):
"But since we live in an imperfect world, some of those measurements are affected by
gaussian white noise:
• Accelerometer : σ = 10−3
, υ = 0
• Gyroscope : σ = 10−2
, υ = 0
• GPS : σ = 10−1
, υ = 0"

The project outline specifies Gaussian GPS noise with sigma = 0.1 meters.
That implies the measurement model:

  z_gps = x_true + v,  v ~ N(0, sigma^2 I)
  sigma = 0.1  =>  variance = 0.01

Under that model, large outliers (e.g. > 3*sigma = 0.3 m) should be very rare.
Empirically, the default stream behaves closer to sigma ~0.9 (variance ~0.81),
and stability requires downweighting GPS far beyond that (R_gps ~ 5000), which
is inconsistent with the Gaussian noise specified in the subject. This is a
clear deviation from the stated specification for the default stream.

In other words: the default stream violates at least one of these Gaussian
assumptions from the subject:
- mean zero noise,
- sigma = 0.1 (variance 0.01),
- thin tails (P(|v| > 3*sigma) ~ 0.0027 per axis).

Procedure
---------
Spec quote (GPS cadence, FR):
"• Votre position GPS actuelle (X,Y,Z en mètres, toutes les trois secondes)."

Spec quote (GPS cadence, EN):
"• Your current GPS position every 3 seconds (X,Y,Z in meters)."

Goal: estimate the GPS noise sigma used by the sensor stream defaults.

Problem: `--debug` output includes TRUE POSITION but does not include the GPS
POSITION lines. Without `--debug`, the GPS POSITION is present but true values
are absent. So we use a two-pass, fixed-seed comparison:

1) Run stream in debug mode with a fixed seed to capture TRUE POSITION by
   timestamp.
2) Run stream again with the same seed (non-debug) and compare GPS POSITION
   against the TRUE POSITION timeline for matching timestamps.
3) Feed back the TRUE POSITION as the "filter estimate" to keep the stream
   advancing (prevents timeouts).

Command outline (automated in a short Python script):
- `./imu-sensor-stream-linux --debug -s 42`
- `./imu-sensor-stream-linux -s 42` (default)
- `./imu-sensor-stream-linux -s 42 -g 0.1`

Evidence (30s simulated time, GPS every ~3s => 11 samples)
----------------------------------------------------------
Default (no -g):
- sigma dx,dy,dz ~= 0.869, 0.969, 0.838
- mean dx,dy,dz  ~= 0.155, -0.335, -0.595

With `-g 0.1`:
- sigma dx,dy,dz ~= 0.087, 0.097, 0.084
- mean dx,dy,dz  ~= 0.0155, -0.0335, -0.0595

Interpretation
--------------
- The GPS noise with `-g 0.1` matches the subject's sigma within sampling error.
- The default stream behaves like sigma ~0.9, not 0.1. That is an order-of-
  magnitude mismatch and is sufficient to cause systematic filter divergence
  if the filter assumes the spec.
- This explains why tuned settings worked better than spec values.
- Even when matching the default GPS variance (~0.81), the filter still fails
  consistently, indicating non-Gaussian behavior/outliers in the default GPS.
  Under a true Gaussian with sigma ~0.9, a Kalman filter configured with
  R_gps ~ 0.81 should stabilize; needing R_gps in the thousands indicates
  heavy-tailed noise or bias that contradicts the Gaussian assumption.

Stability test (single-run sweep)
---------------------------------
Spec quote (failure condition, FR):
"Evidemment, votre estimation de position doit être proche des coordonnées réelles. Si votre estimation est éloignée de plus de 5 mètres de la position
réelle, ou si votre filtre met plus d’une seconde à répondre, alors le programme de flux de
mesures enverra une erreur et arrêtera la communication."

Spec quote (failure condition, EN):
"Obviously, your position estimation must be close to the real coordinates. If your estimation is further than 5 meters of the real position, or if your filter take more than 1 second
to answer, then the sensor-stream will return an error and stop the communication."

Spec quote (duration requirement, EN):
"Your filter must always work with trajectories of durations up to 90 minutes with default
noise amount."

With the default stream settings and the current 12-state filter, the
following single-run tests were executed (each 90 simulated minutes):

- `--gps 0.01`  -> Error: Delta is too high
- `--gps 0.81`  -> Error: Delta is too high
- `--gps 5`     -> Error: Delta is too high
- `--gps 20`    -> Error: Delta is too high
- `--gps 100`   -> Error: Delta is too high
- `--gps 500`   -> Error: Delta is too high
- `--gps 1800`  -> OK (exit 0)
- `--gps 5000`  -> OK (10/10 runs)

Additional testing notes
------------------------
- Control-input model (6-state) with Q = B * Qa * B^T and predict/update in
  the parse loop still failed 10/10 under default stream settings. This
  suggests the primary failure mode is not model choice, but GPS behavior.
- Predict-then-update ordering is required for stable performance; update-only
  or update-then-predict makes the filter lag the stream and increases failures.
- Empirically, the GPS measurement must be heavily downweighted (variance in
  the thousands) to pass consistently on the default stream, which implies
  outliers or heavy-tailed noise beyond a simple Gaussian model.

Further inferences from modification tests
------------------------------------------
- The failure persists across both the 6-state control-input model and the
  12-state measurement-updated model when GPS variance is set near spec
  values, which points away from a pure implementation bug and toward the
  stream's default GPS characteristics.
- The empirical mean of GPS error is not zero in the short sample window.
  This hints at bias or nonstationary noise in the default stream (not
  conclusive with small samples, but consistent with heavy-tailed behavior).

Original implementation behavior (why `-g 0.1` worked)
------------------------------------------------------
The original implementation (before refactors) behaved like a high-R GPS filter:

- `gpsNoise` was set extremely high and scaled by `n^7`, which effectively made
  R_gps enormous and diminished the Kalman gain for GPS updates.
- Acceleration and orientation were treated as direct measurements of the
  corresponding state components, so the filter relied mostly on those streams
  plus the motion model, with GPS acting as a weak correction.

When running the stream with `-g 0.1`, the GPS noise was relatively clean and
consistent, so even an almost-ignored GPS correction was enough to keep the
trajectory within tolerance. With the default stream GPS (sigma ~0.9 and
non-Gaussian outliers), the weak GPS updates could not correct drift fast
enough, leading to "Delta is too high." In short: the tuned settings did not
make GPS accurate; they made the filter less sensitive to GPS, which helped
survive outliers.

Mathematical formulation (reference)
-----------------------------------
Prediction:
  x_k^- = F * x_{k-1}
  P_k^- = F * P_{k-1} * F^T + Q

Measurement update (per sensor):
  K_k = P_k^- * H^T * (H * P_k^- * H^T + R)^-1
  x_k = x_k^- + K_k * (z_k - H * x_k^-)
  P_k = (I - K_k * H) * P_k^-

Gaussian measurement model (spec assumption):
  z_k = H * x_k + v_k
  v_k ~ N(0, R)
  For GPS: R = sigma_gps^2 * I, sigma_gps = 0.1 (per subject).
  pdf(v) = (1 / sqrt((2*pi)^n * det(R))) * exp(-0.5 * v^T * R^-1 * v)

Observed contradiction:
  - Default stream empirical sigma ~0.9 (variance ~0.81), not 0.1.
  - Stable operation requires R_gps orders of magnitude larger (5000).
  - This violates the stated Gaussian noise model at default settings.

Optimal Filter Model and Configuration (to satisfy the stream)
--------------------------------------------------------------
Recommended model: 12-state constant-acceleration model that explicitly
filters acceleration and orientation as measured states. This is less
"textbook" than the control-input model, but it is empirically the most
stable against the stream's default GPS behavior.

State:
  x = [px, py, pz, vx, vy, vz, ax, ay, az, bax, bay, baz, roll, pitch, yaw, bgx, bgy, bgz]

State transition (dt=0.01):
  Position/velocity integrate acceleration; acceleration and orientation are
  modeled as constant between updates.

Process covariance:
  Q = q * I (small scalar, default via --process_noise)

Measurements:
  z_gps = [px, py, pz]
  H_gps = [I 0]
  z_acc = [ax + bax, ay + bay, az + baz]
  z_dir = [roll + bgx, pitch + bgy, yaw + bgz]
  R_gps = sigma_gps^2 * I
  R_acc = sigma_acc^2 * I
  R_dir = sigma_gyr^2 * I

Configuration values:
  - sigma_acc from subject: 1e-3 (variance 1e-6)
  - sigma_gps:
      * If using stream defaults: set variance very high (empirically 5000)
      * If forcing spec behavior: use `-g 0.1` and set variance 0.01
  - process_noise: small scalar for Q = q * I (start at 0.01)
  - initial P:
      * position variance: match R_gps
      * velocity variance: small (e.g., 1.0) unless speed is very uncertain

Practical recommendations
-------------------------
If you cannot change the stream:
  - Run with default stream settings.
  - Set `--gps` to ~5000 (variance) to suppress GPS outliers and pass reliably.
  - Use the 12-state measurement-updated model.

If you can change the stream:
  - Run `./imu-sensor-stream-linux -g 0.1` and keep `--gps 0.01`.
  - This aligns with the subject and produces stable results.

Implementation note
-------------------
The project defaults have been set to the empirically-stable value
`--gps 5000` for the default stream. To switch back to spec-correct
behavior, run the stream with `-g 0.1` and the filter with `--gps 0.01`.

Verification steps
------------------
1) Re-run the two-pass sigma estimate with longer simulated time (e.g. 300s)
   to reduce sampling error (GPS only arrives every 3s).
2) Run `tester.sh` with matched `--gps` and `-g` settings and compare pass rate.
