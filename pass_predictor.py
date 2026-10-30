#!/usr/bin/env python3
"""
Satellite Orbit Ground Track & Overhead Pass Tool
SEDS-TASK-CUB-SRILANKA-012  |  Pod 2C: Comms, Ground Station & Orbital Mechanics
Author: Bismah Farid (IST Chapter)

Standard library only (no external dependencies).

Usage:
    python pass_predictor.py                      # FloripaSat-1 (built-in TLE)
    python pass_predictor.py --tle my_tle.txt     # 3-line TLE file (name + line1 + line2)
    python pass_predictor.py --hours 48 --mask 10 --track track.csv
    python pass_predictor.py --selftest

Model: TLE mean elements (Kozai -> Brouwer conversion) -> Kepler orbit with J2 secular
drift (RAAN, arg. of perigee, mean anomaly) -> ECI -> ECEF (GMST rotation) -> WGS-84
geodetic and topocentric azimuth/elevation.

ACCURACY NOTE: This is a simplified propagator (no drag, no SGP4 periodic terms).
Measure its error against SGP4 with compare_sgp4.py; the +/-0.5 s NORAD target
needs full SGP4.
"""
import argparse
import math
import sys
from datetime import datetime, timedelta, timezone

# ---- Constants (WGS-84 / EGM) ------------------------------------------------
MU = 398600.4418            # km^3/s^2
RE = 6378.137               # km
F = 1.0 / 298.257223563
E2 = F * (2 - F)
J2 = 1.08262668e-3
D2R = math.pi / 180.0
# Orbit-model constants matching the TLE/SGP4 convention (WGS-72)
MU72 = 398600.8             # km^3/s^2
RE72 = 6378.135             # km
J2_72 = 1.082616e-3

# ---- Ground stations (lat deg, lon deg east, altitude km, UTC offset hours) --
STATIONS = {
    "Islamabad": (33.6844, 73.0479, 0.540, 5.0),
    "Colombo":   (6.9271, 79.8612, 0.010, 5.5),
}

# ---- Built-in default satellite: FloripaSat-1 (NORAD 44885) ------------------
FLORIPASAT_TLE = (
    "FLORIPASAT-1",
    "1 44885U 19093G   26302.90796312  .00001457  00000-0  10698-3 0  9995",
    "2 44885  97.7923  38.2356 0010705 340.9320  19.1504 15.05050446372900",
)


# ---- TLE parsing -------------------------------------------------------------
def tle_checksum_ok(line):
    s = 0
    for ch in line[:68]:
        if ch.isdigit():
            s += int(ch)
        elif ch == "-":
            s += 1
    return s % 10 == int(line[68])


class Tle:
    def __init__(self, name, l1, l2):
        if not (l1.startswith("1 ") and l2.startswith("2 ")):
            raise ValueError("Invalid TLE line prefixes")
        if len(l1) < 69 or len(l2) < 69:
            raise ValueError("TLE lines must be 69 characters")
        l1, l2 = l1.rstrip("\r\n"), l2.rstrip("\r\n")
        if not (tle_checksum_ok(l1) and tle_checksum_ok(l2)):
            print("WARNING: TLE checksum mismatch - check the TLE text for typos.",
                  file=sys.stderr)
        self.name = name.strip()
        yr = int(l1[18:20])
        yr += 2000 if yr < 57 else 1900
        day = float(l1[20:32])
        self.epoch = datetime(yr, 1, 1, tzinfo=timezone.utc) + timedelta(days=day - 1)
        self.inc = float(l2[8:16]) * D2R
        self.raan = float(l2[17:25]) * D2R
        self.ecc = float("0." + l2[26:33].strip())
        self.argp = float(l2[34:42]) * D2R
        self.ma = float(l2[43:51]) * D2R
        n_kozai = float(l2[52:63]) * 2 * math.pi / 86400.0   # rad/s, TLE (Kozai) mean motion
        # Convert Kozai -> Brouwer mean motion and semi-major axis (as SGP4 does)
        cosi = math.cos(self.inc)
        x3thm1 = 3 * cosi * cosi - 1
        beta2 = 1 - self.ecc ** 2
        a1 = (MU72 / n_kozai ** 2) ** (1.0 / 3.0) / RE72          # Earth radii
        d1 = 0.75 * J2_72 * x3thm1 / (a1 ** 2 * beta2 ** 1.5)
        a0 = a1 * (1 - d1 / 3 - d1 ** 2 - 134.0 / 81.0 * d1 ** 3)
        d0 = 0.75 * J2_72 * x3thm1 / (a0 ** 2 * beta2 ** 1.5)
        self.n = n_kozai / (1 + d0)                               # Brouwer mean motion, rad/s
        self.a = a0 / (1 - d0) * RE72                             # km
        p = self.a * beta2
        k = J2_72 * (RE72 / p) ** 2
        c = cosi
        self.raan_dot = -1.5 * self.n * k * c
        self.argp_dot = 0.75 * self.n * k * (5 * c * c - 1)
        self.ma_dot = self.n * (1 + 0.75 * k * math.sqrt(beta2) * x3thm1)


def read_tle_file(path):
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    if len(lines) == 2:
        lines.insert(0, "SATELLITE")
    return Tle(lines[0], lines[1], lines[2])


# ---- Time / frames -----------------------------------------------------------
def julian_date(dt):
    dt = dt.astimezone(timezone.utc)
    y, m = dt.year, dt.month
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    day = dt.day + (dt.hour + (dt.minute + (dt.second + dt.microsecond / 1e6) / 60) / 60) / 24
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + day + b - 1524.5


def gmst(dt):
    jd = julian_date(dt)
    t = (jd - 2451545.0) / 36525.0
    g = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * t * t
    return (g % 360.0) * D2R


def solve_kepler(M, e):
    E = M if e < 0.8 else math.pi
    for _ in range(30):
        d = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E -= d
        if abs(d) < 1e-12:
            break
    return E


def propagate_ecef(tle, dt):
    """Return satellite ECEF position (km) at datetime dt."""
    t = (dt - tle.epoch).total_seconds()
    M = (tle.ma + tle.ma_dot * t) % (2 * math.pi)
    raan = tle.raan + tle.raan_dot * t
    argp = tle.argp + tle.argp_dot * t
    e, a, i = tle.ecc, tle.a, tle.inc
    E = solve_kepler(M, e)
    xp = a * (math.cos(E) - e)
    yp = a * math.sqrt(1 - e * e) * math.sin(E)
    cw, sw = math.cos(argp), math.sin(argp)
    cr, sr = math.cos(raan), math.sin(raan)
    ci, si = math.cos(i), math.sin(i)
    x = (cr * cw - sr * sw * ci) * xp + (-cr * sw - sr * cw * ci) * yp
    y = (sr * cw + cr * sw * ci) * xp + (-sr * sw + cr * cw * ci) * yp
    z = (sw * si) * xp + (cw * si) * yp
    g = gmst(dt)
    cg, sg = math.cos(g), math.sin(g)
    return (cg * x + sg * y, -sg * x + cg * y, z)


def geodetic_to_ecef(lat_deg, lon_deg, alt_km):
    lat, lon = lat_deg * D2R, lon_deg * D2R
    n = RE / math.sqrt(1 - E2 * math.sin(lat) ** 2)
    return ((n + alt_km) * math.cos(lat) * math.cos(lon),
            (n + alt_km) * math.cos(lat) * math.sin(lon),
            (n * (1 - E2) + alt_km) * math.sin(lat))


def ecef_to_geodetic(r):
    x, y, z = r
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1 - E2))
    for _ in range(10):
        n = RE / math.sqrt(1 - E2 * math.sin(lat) ** 2)
        alt = p / math.cos(lat) - n
        lat = math.atan2(z, p * (1 - E2 * n / (n + alt)))
    n = RE / math.sqrt(1 - E2 * math.sin(lat) ** 2)
    alt = p / math.cos(lat) - n
    return lat / D2R, ((lon / D2R + 180) % 360) - 180, alt


def look_angles(sat, lat_deg, lon_deg, alt_km):
    """Return (azimuth deg, elevation deg, range km)."""
    gs = geodetic_to_ecef(lat_deg, lon_deg, alt_km)
    d = (sat[0] - gs[0], sat[1] - gs[1], sat[2] - gs[2])
    lat, lon = lat_deg * D2R, lon_deg * D2R
    sl, cl, so, co = math.sin(lat), math.cos(lat), math.sin(lon), math.cos(lon)
    east = -so * d[0] + co * d[1]
    north = -sl * co * d[0] - sl * so * d[1] + cl * d[2]
    up = cl * co * d[0] + cl * so * d[1] + sl * d[2]
    rng = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
    el = math.asin(up / rng) / D2R
    az = math.atan2(east, north) / D2R % 360
    return az, el, rng


# ---- Pass prediction ---------------------------------------------------------
def elevation_at(tle, st, dt):
    lat, lon, alt, _ = st
    return look_angles(propagate_ecef(tle, dt), lat, lon, alt)[1]


def refine_crossing(tle, st, t_lo, t_hi, mask):
    """Bisection to 0.01 s for elevation == mask between t_lo and t_hi."""
    f_lo = elevation_at(tle, st, t_lo) - mask
    for _ in range(40):
        mid = t_lo + (t_hi - t_lo) / 2
        f_mid = elevation_at(tle, st, mid) - mask
        if (f_mid > 0) == (f_lo > 0):
            t_lo, f_lo = mid, f_mid
        else:
            t_hi = mid
        if (t_hi - t_lo).total_seconds() < 0.01:
            break
    return t_lo + (t_hi - t_lo) / 2


def peak_elevation(tle, st, t_a, t_b):
    """Golden-section search for max elevation in [t_a, t_b]."""
    gr = (math.sqrt(5) - 1) / 2
    a, b = 0.0, (t_b - t_a).total_seconds()
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc = elevation_at(tle, st, t_a + timedelta(seconds=c))
    fd = elevation_at(tle, st, t_a + timedelta(seconds=d))
    for _ in range(40):
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = elevation_at(tle, st, t_a + timedelta(seconds=c))
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = elevation_at(tle, st, t_a + timedelta(seconds=d))
    tm = t_a + timedelta(seconds=(a + b) / 2)
    return tm, elevation_at(tle, st, tm)


def find_passes(tle, st, start, hours, mask=5.0, step_s=30):
    end = start + timedelta(hours=hours)
    passes = []
    t = start
    prev = elevation_at(tle, st, t) - mask
    aos = start if prev > 0 else None
    while t < end:
        t2 = t + timedelta(seconds=step_s)
        cur = elevation_at(tle, st, t2) - mask
        if prev <= 0 < cur:
            aos = refine_crossing(tle, st, t, t2, mask)
        elif prev > 0 >= cur and aos is not None:
            los = refine_crossing(tle, st, t, t2, mask)
            tm, emax = peak_elevation(tle, st, aos, los)
            passes.append({
                "aos": aos, "los": los, "tca": tm, "max_el": emax,
                "aos_az": look_angles(propagate_ecef(tle, aos), st[0], st[1], st[2])[0],
                "los_az": look_angles(propagate_ecef(tle, los), st[0], st[1], st[2])[0],
                "duration": (los - aos).total_seconds(),
            })
            aos = None
        t, prev = t2, cur
    return passes


def ground_track(tle, start, hours, step_s=60):
    out, t, end = [], start, start + timedelta(hours=hours)
    while t <= end:
        lat, lon, alt = ecef_to_geodetic(propagate_ecef(tle, t))
        out.append((t, lat, lon, alt))
        t += timedelta(seconds=step_s)
    return out


# ---- Reporting ---------------------------------------------------------------
def fmt(dt, off):
    return (dt + timedelta(hours=off)).strftime("%d-%b %H:%M:%S")


def report(tle, start, hours, mask):
    print(f"Satellite : {tle.name}")
    print(f"TLE epoch : {tle.epoch:%Y-%m-%d %H:%M:%S} UTC")
    print(f"Window    : {start:%Y-%m-%d %H:%M:%S} UTC  +{hours:g} h   (mask {mask:g} deg)")
    print(f"Orbit     : a={tle.a:.1f} km, e={tle.ecc:.5f}, i={tle.inc/D2R:.3f} deg, "
          f"period={2*math.pi/tle.n/60:.2f} min")
    for name, st in STATIONS.items():
        off = st[3]
        label = "PKT" if off == 5.0 else "IST-LK"
        print(f"\n== {name} ({st[0]:.4f}N, {st[1]:.4f}E)  times in UTC+{off:g} ({label}) ==")
        ps = find_passes(tle, st, start, hours, mask)
        if not ps:
            print("  No passes in window.")
            continue
        print(f"  {'AOS':<16}{'TCA':<16}{'LOS':<16}{'MaxEl':>6}{'Dur(s)':>8}{'AOSaz':>7}{'LOSaz':>7}")
        for p in ps:
            print(f"  {fmt(p['aos'],off):<16}{fmt(p['tca'],off):<16}{fmt(p['los'],off):<16}"
                  f"{p['max_el']:>6.1f}{p['duration']:>8.0f}{p['aos_az']:>7.0f}{p['los_az']:>7.0f}")


def selftest():
    # 1. frame sanity: station round trip
    lat, lon, alt = ecef_to_geodetic(geodetic_to_ecef(33.6844, 73.0479, 0.54))
    assert abs(lat - 33.6844) < 1e-6 and abs(lon - 73.0479) < 1e-6 and abs(alt - 0.54) < 1e-6
    # 2. satellite directly overhead => elevation 90
    gs = geodetic_to_ecef(10, 20, 0)
    n = math.sqrt(sum(c * c for c in gs))
    sat = tuple(c / n * (n + 500) for c in gs)
    assert abs(look_angles(sat, 10, 20, 0)[1] - 90) < 0.5
    # 3. GMST at J2000 epoch (2000-01-01 12:00 UTC) ~ 280.46 deg
    assert abs(gmst(datetime(2000, 1, 1, 12, tzinfo=timezone.utc)) / D2R - 280.46) < 0.01
    # 4. TLE + altitude bounds for FloripaSat-1
    t = Tle(*FLORIPASAT_TLE)
    for k in range(0, 5400, 300):
        la, lo, al = ecef_to_geodetic(propagate_ecef(t, t.epoch + timedelta(seconds=k)))
        assert abs(la) <= 90.0 and 500 < al < 700, (la, al)
    print("Self-test passed (frames, overhead geometry, GMST, TLE bounds).")


def main():
    ap = argparse.ArgumentParser(description="Satellite ground track & overhead pass tool")
    ap.add_argument("--tle", help="TLE file (name, line1, line2)")
    ap.add_argument("--hours", type=float, default=24.0)
    ap.add_argument("--mask", type=float, default=5.0, help="elevation mask, degrees")
    ap.add_argument("--start", help="UTC start 'YYYY-MM-DD HH:MM' (default: now)")
    ap.add_argument("--track", help="write ground track CSV to this path")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    tle = read_tle_file(a.tle) if a.tle else Tle(*FLORIPASAT_TLE)
    if a.start:
        start = datetime.strptime(a.start, "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
    else:
        start = datetime.now(timezone.utc)
        if abs((start - tle.epoch).days) > 30:
            print(f"WARNING: TLE is {abs((start - tle.epoch).days)} days old; "
                  "accuracy is poor. Update FLORIPASAT_TLE or use --tle.\n", file=sys.stderr)
    report(tle, start, a.hours, a.mask)
    if a.track:
        with open(a.track, "w") as f:
            f.write("utc,lat_deg,lon_deg,alt_km\n")
            for t, la, lo, al in ground_track(tle, start, a.hours):
                f.write(f"{t:%Y-%m-%dT%H:%M:%SZ},{la:.4f},{lo:.4f},{al:.1f}\n")
        print(f"\nGround track written to {a.track}")


if __name__ == "__main__":
    main()
