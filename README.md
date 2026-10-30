# Satellite Orbit Ground Track & Overhead Pass Tool
SEDS-TASK-CUB-SRILANKA-012 | Bismah Farid, IST Chapter | Satellite: FloripaSat-1 (NORAD 44885)

Run: `python pass_predictor.py` (Python 3, standard library only). Self-test: `--selftest`.

Predicts passes over Islamabad and Colombo and writes the ground track to CSV.
Method: TLE (epoch 2026-10-29) -> Kepler orbit with J2 drift -> ECEF -> azimuth/elevation.
TLE matches the latest published TLE for NORAD 44885.

Limitation: no drag or SGP4 terms, so pass times can differ from official ephemeris by
seconds to tens of seconds; the +/-0.5 s target needs full SGP4 (next step).
