"""Validation on earlier CZ-8A Guowang flights with the same vehicle model.

For each flight with RPHI zones in ../notam/notam_data.json and a published T0, fit the lower stages (kick, core pitch
rate, azimuth) to its two zones, then take the inertial plane at S2 ignition.  Its node relative to GMST(T0), K_L, is a
prediction from NOTAM geometry alone; the measured K comes from the flight's catalog elements (j2_plan.py).  The test
is independent: no TLE enters the zone fit.  Also reported: the core pitch rate the zones need (about 0 means the
gravity turn simply continues, i.e. the masses are consistent with both zones under one kick).

    ~/.venvs/astro313/bin/python validate.py   -> output/validation.json
"""
import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np

import cz8a as C

FLIGHTS = {"Y7": ("2026-007", datetime(2026, 1, 13, 15, 25)), "Y10": ("2026-176", datetime(2026, 8, 4, 8, 52)),
           "Y11": ("2026-221", datetime(2026, 9, 23, 13, 32)), "Y12": ("2026-243", datetime(2026, 10, 25, 2, 3))}


def zones_of(flight):
    d = json.load(open(Path("~/launch/notam/notam_data.json").expanduser()))
    out = []

    def walk(o):
        if isinstance(o, dict):
            if o.get("vehicle") == "CZ-8A" and o.get("type") == "zone" and o.get("flightNum") == flight:
                ring = [tuple(p) for p in o["coords"][:-1]]
                c = C.centroid(ring)
                m = re.search(r"([AB]\d{4}/\d{2})", o["name"])
                out.append(dict(id=m.group(1) if m else o["name"], ring=ring, centroid=c,
                                range_km=float(C.gc_dist(C.PAD, c)), bearing=float(C.bearing(C.PAD, c))))
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(d)
    return sorted(out, key=lambda z: z["range_km"])


def main():
    plan = json.loads((C.OUT / "j2_plan.json").read_text())
    rows = []
    cases = [(k, v[0], v[1]) for k, v in FLIGHTS.items()]
    cases.append(("Y13", "model", datetime(2026, 11, 20, 11, 59)))
    for fl, g, t0 in cases:
        fr = C.Frame(t0)
        z = C.load_zones() if fl == "Y13" else zones_of(fl)
        lw = C.fit_lower(fr, z, x0=(0.5128, 0.05, z[0]["bearing"] - 0.7))
        _, y = C.s2_ignition_state(lw)
        e = C.elements(y[:3], y[3:6])
        f = C.falls(fr, lw)
        k_l = (e["raan"] - np.degrees(fr.g0)) % 360
        k_m = plan["groups"][g]["K_deg"] if g in plan["groups"] else plan["K_y13_deg"]
        rows.append(dict(flight=fl, group=g, t0=t0.isoformat(), zone_bearing=[zz["bearing"] for zz in z],
                         zone_range_km=[zz["range_km"] for zz in z], az=lw["az"], kick=lw["kick"], b_core=lw["b_core"],
                         resid_km=lw["fit_resid_km"], i_lower=e["i"], K_lower=k_l, K_measured=k_m,
                         dK=float(C.WE * 0 + ((k_l - k_m + 180) % 360 - 180)),
                         booster_km=f["boosters"]["range_km"], core_km=f["core"]["range_km"]))
        r = rows[-1]
        print(f"{fl:4s} zones {r['zone_range_km'][0]:.0f}/{r['zone_range_km'][1]:.0f} km brg {r['zone_bearing'][0]:.2f}"
              f" | az {r['az']:.2f} kick {r['kick']:.4f} b {r['b_core']:+.3f} resid {max(map(abs, r['resid_km'])):.1f} km"
              f" | i_L {r['i_lower']:.2f} K_L {r['K_lower']:.3f} K_meas {r['K_measured']:.3f} dK {r['dK']:+.3f}")
    (C.OUT / "validation.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__" and len(__import__("sys").argv) == 1:
    main()


# ---------------------------------------------------------------- phase validation (added 2026-11-21)
PHASE_CASES = {  # flight: (group, T0, measured K, first payload NORAD ids, profile)
    "Y7": ("2026-007", datetime(2026, 1, 13, 15, 25), 310.110, range(67435, 67444), "direct"),
    "Y10": ("2026-176", datetime(2026, 8, 4, 8, 52), 308.470, range(100203, 100212), "longcoast"),
    "Y12": ("2026-243", datetime(2026, 10, 25, 2, 3), 309.696, range(101185, 101194), "direct"),
}


def catalog_first(ids):
    """First catalog TLE of each payload -> (median mean altitude km, median phase lead at T0 helper rows)."""
    H = []
    for f in ("../launch_traj_12/data/hist_joiners.json", "../launch_traj_12/data/hist_2026-007.json",
              "data/hist_2026-243.json"):
        H += json.load(open(C.HERE / f))
    first = {}
    for r in sorted(H, key=lambda r: r["EPOCH"]):
        if int(r["NORAD_CAT_ID"]) in ids:
            first.setdefault(int(r["NORAD_CAT_ID"]), r)
    return list(first.values())


def phase_case(flight, profile=None):
    from sgp4.api import Satrec
    from nrol97 import fit_tle
    g, t0, K, ids, prof = PHASE_CASES[flight]
    prof = profile or prof
    rows = catalog_first(set(ids))
    sats = [Satrec.twoline2rv(r["TLE_LINE1"], r["TLE_LINE2"]) for r in rows]
    lead_cat = float(np.median([C.lead_at_t0(s, t0) for s in sats]))
    alt = float(np.median([(C.MU / (s.no_kozai / 60) ** 2) ** (1 / 3) - C.RE for s in sats]))
    fr = C.Frame(t0)
    z = zones_of(flight)
    lw = C.fit_lower(fr, z, x0=(0.5128, 0.05, z[0]["bearing"] - 0.7))
    nh = C.plane_raan((np.degrees(fr.g0) + K) % 360)
    r_apo = C.RE + alt + 6.0
    for _ in range(3):
        if prof == "direct":
            s2, ab = C.solve_direct(fr, lw, nh, r_apo)
        else:
            s2 = C.solve_s2(fr, lw, nh, C.RE + 200.0, r_apo, x0=[0.5122, -0.1184, 0.1606, -0.0776, 4.3473])
            ab = C.apogee_burn(fr, s2, nh)
        sat = C.march(ab["t_cut"], ab["arc"].Y[-1], ab["t_cut"] + 13000, 5.0, C.M_FREE, np.zeros(22))
        tle, _, _ = fit_tle(fr, sat, ab["t_cut"], ab["t_cut"] + 13000, 99999, "26999A", "X")
        s = Satrec.twoline2rv(tle[-2], tle[-1])
        m_alt = (C.MU / (s.no_kozai / 60) ** 2) ** (1 / 3) - C.RE
        if abs(m_alt - alt) < 0.1:
            break
        r_apo += alt - m_alt
    lead = C.lead_at_t0(s, t0)
    n = np.degrees(s.no_kozai)
    return dict(flight=flight, group=g, profile=prof, mean_alt_km=alt, lead_catalog=lead_cat, lead_model=lead,
                diff_deg=lead - lead_cat, diff_s=(lead - lead_cat) / n * 60.0,
                relight_s=ab["t_ign"], s2_left_kg=C.VEH.s2_prop - s2["prop_used"] - ab["prop"])


def phase_main():
    rows = []
    for fl, prof in (("Y10", "longcoast"), ("Y10", "direct"), ("Y12", "direct"), ("Y12", "longcoast"),
                     ("Y7", "direct")):
        try:
            r = phase_case(fl, prof)
        except RuntimeError as err:      # direct sequence without a solution for this flight's zones and load
            rows.append(dict(flight=fl, profile=prof, error=str(err)))
            print(f"{fl:4s} {prof:9s} {err}", flush=True)
            continue
        rows.append(r)
        print(f"{fl:4s} {prof:9s} park {r['mean_alt_km']:7.1f} km | lead catalog {r['lead_catalog']:+7.2f} model "
              f"{r['lead_model']:+7.2f} -> {r['diff_deg']:+6.2f} deg ({r['diff_s']:+6.1f} s) | relight T+"
              f"{r['relight_s']:.0f} s, S2 left {r['s2_left_kg']:.0f} kg", flush=True)
    (C.OUT / "phase_validation.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__" and len(__import__("sys").argv) > 1 and __import__("sys").argv[1] == "phase":
    phase_main()
