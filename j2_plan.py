"""T0 and target plane of CZ-8A Y13 from the Guowang plane-joining plan (J2 differential drift).

Site constant K = RAAN_inj(T0) - GMST(T0) and tuning offset D = RAAN_inj(T0) - RAAN_target(T0) are measured on the
recent CZ-8A groups with published T0; Y13 is assumed to join the plane of 2025-295 (CZ-12 Y4), as the second group.
The functions of ../launch_traj_12/j2_tuning.py are imported, not copied.  Also tabulates the NOTAM-to-lift-off lead of
every CZ-8A flight in ../notam/notam_data.json.

    ~/.venvs/astro313/bin/python j2_plan.py        -> output/j2_plan.json
"""
import json
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from astropy.time import Time
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path("~/launch/launch_traj_12").expanduser()))
import j2_tuning as J  # noqa: E402

MU, RE = J.MU, J.RE
# published T0 (UTC): Xinhua / CCTV, GCAT
T0_PUB = {"2026-007": datetime(2026, 1, 13, 15, 25), "2026-176": datetime(2026, 8, 4, 8, 52),
          "2026-221": datetime(2026, 9, 23, 13, 32), "2026-243": datetime(2026, 10, 25, 2, 3)}
FLIGHT = {"2026-007": "Y7", "2026-176": "Y10", "2026-221": "Y11", "2026-243": "Y12"}
TARGET = "2025-295"
DAY = datetime(2026, 11, 20, 11, 0)
T0_NOTAM = datetime(2026, 11, 20, 11, 59)        # zone-1 B time 11:51 + 8 min (pattern below)


def gmst(t):
    return Time(t, scale="utc").sidereal_time("mean", "greenwich").deg


def parse(rows, ids):
    out = []
    for r in rows:
        if int(r["NORAD_CAT_ID"]) not in ids:
            continue
        n = float(r["MEAN_MOTION"]) * 2 * np.pi / 86400.0
        out.append((Time(datetime.fromisoformat(r["EPOCH"])).mjd, float(r["RA_OF_ASC_NODE"]),
                    (MU / n ** 2) ** (1 / 3) - RE, float(r["INCLINATION"])))
    a = np.array(sorted(out))
    return dict(t=a[:, 0], raan=a[:, 1], alt=a[:, 2], inc=a[:, 3])


def target_raan(t):
    """Mean plane of the 2025-295 payloads, propagated from current gp with the measured working-orbit rate."""
    rows = [r for r in json.load(open(HERE / "data" / "gp_2025-295.json")) if r["OBJECT_TYPE"] == "PAYLOAD"]
    vals = [float(r["RA_OF_ASC_NODE"]) + W_RATE * (Time(t).mjd - Time(datetime.fromisoformat(r["EPOCH"])).mjd)
            for r in rows]
    return float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(vals))))) % 360)


def notam_leads():
    d = json.load(open(Path("~/launch/notam/notam_data.json").expanduser()))
    out = {}

    def walk(o):
        if isinstance(o, dict):
            if o.get("vehicle") == "CZ-8A" and o.get("type") == "zone":
                m = re.search(r"B\)\s*(\d{10})\s*C\)\s*(\d{10})", o["description"])
                if m:
                    out.setdefault(o["flightNum"], []).append(m.groups())
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(d)
    return out


hist, _ = J.load_groups()
planes = J.current_planes()
W_RATE, _ = J.working_rate(hist)
hist["2026-243"] = parse(json.load(open(HERE / "data" / "hist_2026-243.json")), set(range(101185, 101194)))
PAIRS = {"2026-176": "2025-187", "2026-221": "2025-285", "2026-243": "2025-311", "2026-007": None}


def main():
    res = {"working_rate_deg_d": W_RATE, "groups": {}}
    for g, t0 in T0_PUB.items():
        h, m0 = hist[g], Time(t0).mjd
        park, t_raise, t_arr, _, _ = J.parking_phase(h, m0)
        r0, rate, rms, n = J.injection_raan(h, m0, park)
        tg = PAIRS[g]
        D = float(J.wrap180(r0 - J.target_raan(tg, m0, hist, planes, W_RATE))) if tg else None
        res["groups"][g] = dict(flight=FLIGHT[g], t0=t0.isoformat(), park_alt_km=float(park), raise_start_d=float(t_raise),
                                arrival_d=None if np.isnan(t_arr) else float(t_arr), raan_inj=float(r0),
                                rate_deg_d=float(rate), K_deg=float((r0 - gmst(t0)) % 360), target=tg, D_deg=D)
    K12, D12 = res["groups"]["2026-243"]["K_deg"], res["groups"]["2026-243"]["D_deg"]
    K1011 = np.mean([res["groups"][g]["K_deg"] for g in ("2026-176", "2026-221")])
    D1011 = np.mean([res["groups"][g]["D_deg"] for g in ("2026-176", "2026-221")])
    sols = {}
    for lab, K, D in (("Y12 plan (park 1050 km)", K12, D12), ("Y10/Y11 plan (park 1100 km)", K1011, D1011)):
        f = lambda s: J.wrap180(gmst(DAY + timedelta(seconds=s)) + K - target_raan(DAY + timedelta(seconds=s)) - D)  # noqa
        s = brentq(f, 0, 7200)
        sols[lab] = dict(K=float(K), D=float(D), T0=(DAY + timedelta(seconds=s)).isoformat(timespec="seconds"))
    leads = {}
    for fl, zs in notam_leads().items():
        g = next((k for k, v in FLIGHT.items() if v == fl), None)
        leads[fl] = dict(B=sorted(b for b, _ in zs), C=sorted(c for _, c in zs))
    res.update(target=TARGET, target_raan_at_notam_t0=target_raan(T0_NOTAM), solutions=sols, notam=leads,
               t0_notam=T0_NOTAM.isoformat(), K_y13_deg=K12,
               raan_y13_at_notam_t0=float((gmst(T0_NOTAM) + K12) % 360),
               D_y13_at_notam_t0=float(J.wrap180(gmst(T0_NOTAM) + K12 - target_raan(T0_NOTAM))))
    out = HERE / "missions" / "Y13_20261120" / "output"   # = cz8a.OUT
    out.mkdir(exist_ok=True)
    (out / "j2_plan.json").write_text(json.dumps(res, indent=1))
    for g, r in res["groups"].items():
        print(f"{r['flight']:4s} {g} T0 {r['t0']} park {r['park_alt_km']:.1f} km K {r['K_deg']:.3f} "
              f"target {r['target']} D {r['D_deg'] if r['D_deg'] is None else round(r['D_deg'], 3)}")
    for k, v in sols.items():
        print(f"{k:30s} T0 {v['T0']}")
    print(f"NOTAM T0 {T0_NOTAM}: 2025-295 RAAN {res['target_raan_at_notam_t0']:.3f}, Y13 RAAN "
          f"{res['raan_y13_at_notam_t0']:.3f}, D {res['D_y13_at_notam_t0']:+.3f}")


if __name__ == "__main__":
    main()
