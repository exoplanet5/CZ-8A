"""Compare the model orbit with the outside element set outside-ref/outside.tle.

That element set is Y12's first catalog TLE of 2026-243E (101189, epoch 2026-10-25 11:59:22.887) moved to Y13:
epoch +26 d 09:57:00, RAAN + the GMST change over that interval, mean anomaly +1.5523 deg (27.5 s).  It reported a
successful catch of the Y13 payload near T+100 min, so it stands for "Y13 flew exactly like Y12".

For each profile the script gives the radial / along-track / cross-track offset of the model from the outside orbit,
the argument-of-latitude and node differences, and the T0 shift each would imply:
    a pure T0 error dt moves the node by +0.2507 deg/min * dt and the phase by -n * dt; the two must agree.

    ~/.venvs/astro313/bin/python compare_outside.py   -> output/outside_comparison.json
"""
import json
from datetime import datetime, timedelta

import numpy as np
from sgp4.api import Satrec

import build_traj as B
import cz8a as C
from nrol97 import jd_of

T0 = B.T0_NOMINAL
WE_DEG_MIN = np.degrees(C.WE) * 60.0


def outside():
    L = (C.MISSION / "outside-ref" / "outside.tle").read_text().splitlines()
    return Satrec.twoline2rv(L[1], L[2]), L


def y12_source():
    rows = sorted([r for r in json.load(open(C.HERE / "data" / "hist_2026-243.json"))
                   if int(r["NORAD_CAT_ID"]) == 101189], key=lambda r: r["EPOCH"])
    return rows[0]


def offsets(frame, model, sat, mins):
    out = []
    n_deg_min = np.degrees(sat.no_kozai)          # rad/min -> deg/min
    for m in mins:
        t = m * 60.0
        _, ro, vo = sat.sgp4(jd_of(frame.utc(t)), 0.0)
        ro, vo = np.array(ro), np.array(vo)
        s = model(t)
        rm, vm = s[:3], s[3:6]
        R = ro / np.linalg.norm(ro)
        W = np.cross(ro, vo)
        W /= np.linalg.norm(W)
        S = np.cross(W, R)
        d = rm - ro
        eo, em = C.elements(ro, vo), C.elements(rm, vm)
        du = float((em["u"] - eo["u"] + 180) % 360 - 180)
        dO = float((em["raan"] - eo["raan"] + 180) % 360 - 180)
        out.append(dict(met_min=m, utc=frame.utc(t).isoformat(timespec="seconds"), dist_km=float(np.linalg.norm(d)),
                        radial_km=float(d @ R), along_km=float(d @ S), cross_km=float(d @ W),
                        along_s=float(d @ S / np.linalg.norm(vo)), du_deg=du, draan_deg=dO,
                        di_deg=float(em["i"] - eo["i"]), dr_km=float(np.linalg.norm(rm) - np.linalg.norm(ro)),
                        t0_shift_from_phase_s=float(du / n_deg_min * 60.0),
                        t0_shift_from_node_s=float(-dO / WE_DEG_MIN * 60.0)))
    return out


def main():
    sat, lines = outside()
    src = y12_source()
    ep_out = datetime(2026, 1, 1) + timedelta(days=float(lines[1][20:32]) - 1)
    ep_src = datetime.fromisoformat(src["EPOCH"])
    l2, s2 = lines[2].split(), src["TLE_LINE2"].split()
    origin = dict(source="2026-243E (101189) first catalog TLE", source_epoch=src["EPOCH"],
                  epoch_shift=str(ep_out - ep_src), same_i_e_argp_n=(l2[2], l2[4], l2[5], l2[7][:11]) ==
                  (s2[2], s2[4], s2[5], s2[7][:11]), d_mean_anomaly_deg=float(l2[6]) - float(s2[6]))
    plan = json.loads((C.OUT / "j2_plan.json").read_text())
    res = dict(outside_tle=lines, origin=origin, t0=T0.isoformat(),
               lead_at_t0_deg=dict(outside=C.lead_at_t0(sat, T0)), profiles={})
    mins = [40, 70, 80, 90, 100, 110, 120, 150, 200, 240]
    for prof in ("attached", "direct", "longcoast"):
        try:
            sol = B.solve(T0, B.VEH_NOM, plan["K_y13_deg"], prof, verbose=False)
        except RuntimeError as err:      # the direct sequence does not close from the revision 5 lower stages
            res["profiles"][prof] = dict(error=str(err))
            print(f"{prof}: {err}")
            continue
        segs = B.segments(sol)
        arcs = [a for _, a in segs] + [sol["sat"]]

        def model(t, arcs=arcs):
            for a in arcs:
                if a.T[0] <= t <= a.T[-1]:
                    return a(t)
            raise ValueError(t)
        rows = offsets(sol["frame"], model, sat, mins)
        res["profiles"][prof] = dict(offsets=rows, lead_at_t0_deg=sol["lead"], mean_alt_km=sol["mean_alt"],
                                     tle=sol["tle"])
        res["lead_at_t0_deg"][prof] = sol["lead"]
        print(f"{prof}: lead at T0 {sol['lead']:+.2f} deg (outside {res['lead_at_t0_deg']['outside']:+.2f})")
        print("   MET   along km    s   du deg  dRAAN   di    dr km | T0 shift from phase / from node (s)")
        for r in rows:
            print(f"  {r['met_min']:4d} {r['along_km']:8.1f} {r['along_s']:6.1f} {r['du_deg']:+7.3f} "
                  f"{r['draan_deg']:+.3f} {r['di_deg']:+.3f} {r['dr_km']:+6.2f} | "
                  f"{r['t0_shift_from_phase_s']:+7.1f} / {r['t0_shift_from_node_s']:+6.1f}")
    (C.OUT / "outside_comparison.json").write_text(json.dumps(res, indent=1))
    print("origin:", origin)


if __name__ == "__main__":
    main()
