"""
Validation only (NOT part of the deliverable): measure the error of pass_predictor.py
against the SGP4 reference model for the same TLE.
Setup (once):  pip install sgp4
Run:           python compare_sgp4.py
"""
import math
from datetime import datetime, timezone
from sgp4.api import Satrec, jday
import pass_predictor as pp

L1, L2 = pp.FLORIPASAT_TLE[1], pp.FLORIPASAT_TLE[2]
sat = Satrec.twoline2rv(L1, L2)
mine = pp.propagate_ecef


def sgp4_ecef(tle, dt):
    dt = dt.astimezone(timezone.utc)
    jd, fr = jday(dt.year, dt.month, dt.day, dt.hour, dt.minute,
                  dt.second + dt.microsecond / 1e6)
    err, r, v = sat.sgp4(jd, fr)
    if err:
        raise RuntimeError(f"SGP4 error code {err}")
    g = pp.gmst(dt)                       # TEME -> ECEF (z rotation by GMST)
    cg, sg = math.cos(g), math.sin(g)
    return (cg * r[0] + sg * r[1], -sg * r[0] + cg * r[1], r[2])


tle = pp.Tle(*pp.FLORIPASAT_TLE)
start = datetime(2026, 10, 30, 19, 0, tzinfo=timezone.utc)
HOURS = 24

# 1) position error vs time since TLE epoch
print("Position difference, my model vs SGP4:")
for hrs in (0, 6, 12, 24, 36):
    t = tle.epoch + pp.timedelta(hours=hrs)
    a, b = mine(tle, t), sgp4_ecef(tle, t)
    d = math.dist(a, b)
    print(f"  epoch +{hrs:>2} h : {d:7.2f} km  (~{d / 7.6:5.1f} s along-track)")

# 2) pass-time error
for name, st in pp.STATIONS.items():
    pp.propagate_ecef = mine
    p_mine = pp.find_passes(tle, st, start, HOURS)
    pp.propagate_ecef = sgp4_ecef
    p_ref = pp.find_passes(tle, st, start, HOURS)
    pp.propagate_ecef = mine
    print(f"\n{name}: pass-time difference (mine - SGP4)")
    for a in p_mine:
        b = min(p_ref, key=lambda x: abs((x['aos'] - a['aos']).total_seconds()), default=None)
        if b is None:
            continue
        print(f"  AOS {a['aos']:%d-%b %H:%M:%S}Z  dAOS={(a['aos']-b['aos']).total_seconds():+6.1f} s"
              f"  dLOS={(a['los']-b['los']).total_seconds():+6.1f} s"
              f"  dMaxEl={a['max_el']-b['max_el']:+5.2f} deg")

