"""Build the CZ-8A Y13 nominal trajectory: zone fit, S2 two-burn to the 1050 km parking orbit, disposal, outputs.

    ~/.venvs/astro313/bin/python build_traj.py [--t0 2026-11-20T11:59:59] [--payload 6300] [--profile direct]
                                               [--shift S] [--u-apo DEG] [--full-thrust] [--sens]

All profiles share the lower stages of revision 5: fitted to both NOTAM zones and the four ground sightings
(sightings.py), core engines at 77.5 % from T+10 s to booster separation (VEH_NOM); --full-thrust flies the published
371 t vehicle instead (the sightings then miss by up to 2.5 deg).
--profile attached (default): lofted S2 burn 1, a 40.7 min coast to an apogee on the inserted orbit
(outside-ref/outside.tle), burn 2 onto its velocity.  The vehicle runs --shift seconds ahead of that element set; by
default the shift is 0, or the smallest one the S2 load can pay for (S2_RESERVE left after burn 2).
--profile direct: the published CZ-8A sequence, S2 cutoff T+709.4 s on a lofted sub-orbital arc, relight
T+1761.5 s at the apex (Y12's catalog phase supports it).
--profile longcoast: the Y10/Y11 type, perigee cutoff, half-orbit coast, apogee burn near T+59 min.

Writes output/CZ-8A_Y13_{traj.csv, orbit.csv, falls.csv, solution.json, prelaunch.tle, traj.kml}.
The target plane comes from j2_plan.json (j2_plan.py): RAAN = GMST(T0) + K, K from Y12 (2026-243).
"""
import argparse
import csv
import json
from dataclasses import replace
from datetime import datetime, timedelta

import numpy as np
from sgp4.api import Satrec

import cz8a as C
import sightings as SG
from nrol97 import J2, fit_tle

OUT = C.OUT
OUT.mkdir(exist_ok=True)
RP_STAGE = C.RE + 904.6      # Y12 stage perigee after its disposal burn (2026-243K)
T_DEPLOY = 120.0             # s after the apogee burn, ASSUMED
T_DISPOSAL = 600.0           # long-coast profile only: s after deployment, ASSUMED


def y12_reference():
    """Y12 (2026-243) catalog targets: median mean altitude of the nine payloads' first TLEs, and the argument of
    latitude of the stage's apogee (its last burn) at T0, from 2026-243K's first TLE with the J2 apsidal rate."""
    rows = sorted(json.load(open(C.HERE / "data" / "hist_2026-243.json")), key=lambda r: r["EPOCH"])
    first = {}
    for r in rows:
        first.setdefault(int(r["NORAD_CAT_ID"]), r)
    alts = [(C.MU / (float(first[k]["MEAN_MOTION"]) * 2 * np.pi / 86400) ** 2) ** (1 / 3) - C.RE
            for k in range(101185, 101194)]
    rb = first[101194]
    n = float(rb["MEAN_MOTION"]) * 2 * np.pi / 86400
    a, e, i = (C.MU / n ** 2) ** (1 / 3), float(rb["ECCENTRICITY"]), np.radians(float(rb["INCLINATION"]))
    wdot = np.degrees(0.75 * n * J2 * (C.RE / (a * (1 - e * e))) ** 2 * (5 * np.cos(i) ** 2 - 1))     # deg/s
    dt = (datetime.fromisoformat(rb["EPOCH"]) - datetime(2026, 10, 25, 2, 3)).total_seconds()
    argp0 = float(rb["ARG_OF_PERICENTER"]) - wdot * dt
    return dict(mean_alt_km=float(np.median(alts)), stage_apogee_u_deg=float((argp0 + 180) % 360))


Y12 = y12_reference()


T0_NOMINAL = datetime(2026, 11, 20, 11, 59, 59)   # inside the published minute.  Since revision 5 the phase no longer
                                                  # pins it: zero shift closes from ~11:59:41 on (coast_tradeoff.py)
U_APO_ATTACHED = 323.0      # deg: apogee point on the inserted orbit; gives a 40.7 min coast (the user's range is 40-48)
AIM_KM = (0.0, 0.0)         # booster / core aim points along-track from the zone centroids (30 km for the boosters in
                            # revision 4, to lift the fairing over 100 km; the sightings now put it at 135 km)
# Revision 5 nominal vehicle: the core throttled to 77.5 % (the CZ-8 throttle level) from T+10 s to booster separation.
# Only with this much less energy do the four sightings and both zones fit together; the propellant it saves is not
# loaded, so lift-off is 341.3 t against the published 371 t (full thrust: 371.0 t with the 30.2 t S2).
VEH_NOM = replace(C.VEH, core_thr=0.775, t_thr=10.0)
X0_LOWER = (0.3954, 0.5131, 139.236, 135.75, -1.063)   # kick, b_core, az, t_q, b_boost (fit_lower, sight=True)
S2_RESERVE = 130.0          # kg left after burn 2 in the attached profile (disposal uses ~30 kg)
X0_ATTACHED = [5.3, -1.14, 0.043, -0.018, 5.87]


def attach(frame, lw, u_apo, veh, shift=None, reserve=S2_RESERVE, hi=240.0):
    """Attached long coast from the lower stages lw.  shift: seconds the vehicle runs ahead of the inserted element
    set; None solves the smallest shift >= 0 that leaves `reserve` kg of S2 propellant after burn 2.
    -> (s2, ab, t_apo, shift, ins)"""
    from scipy.optimize import brentq
    st = {"x0": X0_ATTACHED}

    def run(sh):
        ins = C.inserted_orbit(shift=sh)
        s2, ab, t_apo = C.solve_attached(frame, lw, ins, u_apo, veh, x0=st["x0"])
        if s2["cost"] > 1e-6:
            raise RuntimeError(f"attached profile: no solution at shift {sh:.1f} s (cost {s2['cost']:.1e})")
        st["x0"] = s2["x"]
        return s2, ab, t_apo, float(sh), ins

    def margin(sh):
        s2, ab = run(sh)[:2]
        return veh.s2_prop - s2["prop_used"] - ab["prop"] - reserve

    if shift is None:
        # walk down from a shift the start values solve, 25 s at a time (continuation), to the sign change
        up = 100.0
        while margin(up) < 0:
            up += 25.0
            if up > hi:
                raise RuntimeError(f"attached profile: the S2 load does not close below a {hi:.0f} s shift")
        lo = up
        while lo > 0 and margin(lo := max(lo - 25.0, 0.0)) >= 0:
            up = lo
        shift = 0.0 if up == 0.0 else brentq(margin, lo, up, xtol=0.05)
    return run(shift)


def solve(t0, veh, k_site, profile="attached", r_cut_alt=200.0, t_ign2=C.T_IGN2_DIRECT, u_apo=U_APO_ATTACHED,
          shift=None, aim=AIM_KM, verbose=True):
    """All profiles fly the same lower stages, fitted to both zones (aimed at the centroids + aim) and the ground
    sightings (fit_lower, sight=True).
    profile 'attached' (nominal, revision 5): lofted burn 1 to a high-perigee transfer, a 40.7 min coast to an apogee
    placed on the inserted orbit (outside-ref/outside.tle) at argument of latitude u_apo, burn 2 ending on that
    orbit's velocity.  The vehicle is `shift` seconds ahead of the element set (None: the smallest shift the S2 load
    pays for, see attach).
    'direct': Y12-type direct ascent (published sequence).  'longcoast': Y10-type perigee cutoff at r_cut_alt,
    half-orbit coast, apogee burn into the Y12 mean altitude."""
    frame = C.Frame(t0)
    zones = C.load_zones()
    lw = C.fit_lower(frame, zones, veh, x0=X0_LOWER, aim=aim, sight=True)
    raan = (np.degrees(frame.g0) + k_site) % 360
    nh = C.plane_raan(raan)
    target = Y12["mean_alt_km"]
    r_apo = C.RE + target + 6.0
    x0 = [0.5122, -0.1184, 0.1606, -0.0776, 4.3473]
    for it in range(5):
        if profile == "attached":
            s2, ab, t_apo_ins, shift, ins = attach(frame, lw, u_apo, veh, shift)
            el_ins = C.elements(*C.sgp4_teme(ins, frame, 0.0))
            raan = el_ins["raan"]
            k_site = (raan - np.degrees(frame.g0)) % 360
            nh = el_ins["hhat"]
        elif profile == "direct":
            s2, ab = C.solve_direct(frame, lw, nh, r_apo, veh, t_ign2=t_ign2)
        else:
            s2 = C.solve_s2(frame, lw, nh, C.RE + r_cut_alt, r_apo, veh, x0=x0)
            x0 = s2["x"]
            ab = C.apogee_burn(frame, s2, nh, veh)
        y = ab["arc"].Y[-1].copy()
        t_dep = ab["t_cut"] + T_DEPLOY
        c = C.march(ab["t_cut"], y, t_dep, 2.0, C.M_FREE, np.zeros(22))
        sat = C.march(t_dep, c.Y[-1], t_dep + 2 * 6500.0, 5.0, C.M_FREE, np.zeros(22))
        tle, rms, mx = fit_tle(frame, sat, t_dep, t_dep + 2 * 6500.0, 99913, "26999A", "CZ-8A Y13 PAYLOAD")
        s = Satrec.twoline2rv(tle[-2], tle[-1])
        mean_alt = (C.MU / (s.no_kozai / 60.0) ** 2) ** (1 / 3) - C.RE
        if verbose:
            print(f"  iter {it}: r_apo-RE {r_apo - C.RE:.2f} km -> mean alt {mean_alt:.2f} km "
                  f"(target {target:.2f}, TLE fit rms {rms:.3f} km)")
        if profile == "attached" or abs(mean_alt - target) < 0.05:
            break
        r_apo += target - mean_alt
    # deployment; stage disposal (retrograde, perigee to the Y12 stage value)
    y_dep = c.Y[-1].copy()
    y_dep[6] -= veh.payload
    if profile == "attached":
        t_dsp = t_dep + 300.0           # disposal by propellant dump soon after deployment, ASSUMED
    elif profile == "direct":
        # at the argument of latitude where Y12's stage made its last burn (its apogee)
        probe = C.march(t_dep, y_dep, t_dep + 6500.0, 2.0, C.M_FREE, np.zeros(22))
        us = np.array([C.elements(Y[:3], Y[3:6])["u"] for Y in probe.Y])
        d = (us - Y12["stage_apogee_u_deg"] + 180) % 360 - 180
        k = int(np.flatnonzero((d[:-1] < 0) & (d[1:] >= 0))[0])
        t_dsp = probe.T[k] - d[k] / (d[k + 1] - d[k]) * (probe.T[k + 1] - probe.T[k])
    else:
        t_dsp = t_dep + T_DISPOSAL
    stage = C.march(t_dep, y_dep, t_dsp, 2.0, C.M_FREE, np.zeros(22))
    dsp = C.disposal(frame, stage.Y[-1], t_dsp, nh, RP_STAGE, veh)
    rb = C.march(dsp["arc"].T[-1], dsp["arc"].Y[-1], dsp["arc"].T[-1] + 2 * 6500.0, 5.0, C.M_FREE, np.zeros(22))
    tle_rb, rms_rb, _ = fit_tle(frame, rb, rb.T[0], rb.T[-1], 99914, "26999K", "CZ-8A Y13 R/B")
    fl = C.falls(frame, lw, veh, s2)
    return dict(frame=frame, zones=zones, lw=lw, nh=nh, raan=raan, s2=s2, ab=ab, coast_dep=c, sat=sat, tle=tle,
                tle_rms=rms, mean_alt=mean_alt, r_apo=r_apo, stage=stage, dsp=dsp, rb=rb, tle_rb=tle_rb, falls=fl,
                t_dep=t_dep, veh=veh, k_site=k_site, r_cut_alt=r_cut_alt, profile=profile, t_ign2=t_ign2,
                lead=C.lead_at_t0(s, t0), u_apo=u_apo if profile == "attached" else None,
                shift=shift if profile == "attached" else None)


def segments(sol):
    """(phase, Arc) in time order for the vehicle carrying the payload, then the stage after deployment."""
    v, lw, s2, ab = sol["veh"], sol["lw"], sol["s2"], sol["ab"]
    a = lw["arc"]
    i_b = int(np.searchsorted(a.T, v.t_bsep))
    seg = [("boost", C.Arc(a.T[:i_b + 1], a.Y[:i_b + 1])), ("core", C.Arc(a.T[i_b:], a.Y[i_b:]))]
    s2a = s2["arc"]
    j = int(np.searchsorted(s2a.T, v.t_s2))
    seg += [("staging", C.Arc(s2a.T[:j + 1], s2a.Y[:j + 1])), ("S2 burn 1", C.Arc(s2a.T[j:], s2a.Y[j:])),
            ("coast", ab["coast"]), ("S2 burn 2", ab["arc"]), ("orbit", sol["coast_dep"])]
    return seg


def sample(arc, t0, t1, step):
    ts = np.unique(np.r_[np.arange(np.ceil(t0 / step) * step, t1, step), t0, t1])
    return [(t, arc(t)) for t in ts if t0 <= t <= t1]


def row(frame, phase, t, s, pad=C.PAD):
    r, v = s[:3], s[3:6]
    lat, lon, h = frame.lla(r, t)
    vr = C.v_rel(r, v)
    rn = np.linalg.norm(r)
    rr = frame.to_ecef(r, t)
    return dict(t=f"{t:.3f}", utc=frame.utc(t).isoformat(timespec="milliseconds"), phase=phase, lat=f"{lat:.6f}",
                lon=f"{lon:.6f}", alt_km=f"{h:.4f}", x=f"{rr[0]:.4f}", y=f"{rr[1]:.4f}", z=f"{rr[2]:.4f}",
                v_inertial=f"{np.linalg.norm(v):.6f}", v_rel=f"{np.linalg.norm(vr):.6f}",
                fpa=f"{np.degrees(np.arcsin(np.clip((r @ vr) / (rn * max(np.linalg.norm(vr), 1e-9)), -1, 1))):.6f}",
                downrange_km=f"{C.gc_dist(pad, (lat, lon)):.4f}")


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def kml(sol, path):
    fr = sol["frame"]

    def coords(arc, t0, t1, step):
        out = []
        for t, s in sample(arc, t0, t1, step):
            la, lo, h = fr.lla(s[:3], t)
            out.append(f"{lo:.5f},{la:.5f},{max(h, 0) * 1e3:.0f}")
        return " ".join(out)

    def line(name, color, c, alt=True):
        mode = "absolute" if alt else "clampToGround"
        return (f"<Placemark><name>{name}</name><Style><LineStyle><color>{color}</color><width>2</width></LineStyle>"
                f"</Style><LineString><altitudeMode>{mode}</altitudeMode><coordinates>{c}</coordinates></LineString>"
                f"</Placemark>")

    parts = []
    for ph, arc in segments(sol):
        parts.append(line(ph, "ff4fc3f7" if "burn" in ph or ph in ("boost", "core") else "ffffffff",
                          coords(arc, arc.T[0], arc.T[-1], 2.0 if arc.T[-1] < 800 else 20.0)))
    parts.append(line("payload orbit (2 rev)", "ff9e9e9e", coords(sol["sat"], sol["sat"].T[0], sol["sat"].T[-1], 30.0)))
    for k, f in sol["falls"].items():
        parts.append(line(f"{k} fall", "ff00a5ff", coords(f["arc"], f["arc"].T[0], f["arc"].T[-1], 2.0)))
        parts.append(f"<Placemark><name>{k} impact T+{f['t']:.0f} s</name><Point><coordinates>{f['lon']:.4f},"
                     f"{f['lat']:.4f},0</coordinates></Point></Placemark>")
    for z in sol["zones"]:
        ring = z["ring"] + [z["ring"][0]]
        c = " ".join(f"{lo:.5f},{la:.5f},0" for la, lo in ring)
        parts.append(f"<Placemark><name>{z['id']} {z['t0']}-{z['t1']}</name><Style><LineStyle><color>ff3d5afe"
                     f"</color><width>2</width></LineStyle><PolyStyle><color>403d5afe</color></PolyStyle></Style>"
                     f"<Polygon><outerBoundaryIs><LinearRing><coordinates>{c}</coordinates></LinearRing>"
                     f"</outerBoundaryIs></Polygon></Placemark>")
    path.write_text('<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document>'
                    f"<name>CZ-8A Y13 {fr.t0:%Y-%m-%d %H:%M} UTC</name>{''.join(parts)}</Document></kml>")


def summary(sol):
    fr, v, lw, s2, ab, dsp = sol["frame"], sol["veh"], sol["lw"], sol["s2"], sol["ab"], sol["dsp"]
    q, mach = C.q_profile(lw["arc"])
    iq = int(np.argmax(q))
    e_tr = C.elements(s2["arc"].r_end, s2["arc"].v_end)
    e_park = C.elements(ab["arc"].r_end, ab["arc"].v_end)
    e_rb = C.elements(dsp["arc"].r_end, dsp["arc"].v_end)
    c_ign, y_ign = C.s2_ignition_state(lw, v)
    ev = {
        "booster_sep": C.state_summary(fr, v.t_bsep, lw["y_bsep"]),
        "core_cutoff": C.state_summary(fr, v.t_ccut, lw["arc"].Y[-1]),
        "s2_ignition": C.state_summary(fr, v.t_s2, y_ign),
        "fairing": C.state_summary(fr, v.t_fair, s2["arc"](v.t_fair)),
        "s2_cutoff_1": C.state_summary(fr, s2["t_cut"], s2["arc"].Y[-1]),
        "s2_ignition_2": C.state_summary(fr, ab["t_ign"], ab["arc"].Y[0]),
        "s2_cutoff_2": C.state_summary(fr, ab["t_cut"], ab["arc"].Y[-1]),
        "deployment": C.state_summary(fr, sol["t_dep"], sol["coast_dep"].Y[-1]),
        "disposal": C.state_summary(fr, dsp["t_ign"], dsp["arc"].Y[0]),
    }
    for k in ev:
        ev[k]["u_deg"] = None
    for k, (t, s) in {"s2_cutoff_1": (s2["t_cut"], s2["arc"].Y[-1]), "s2_cutoff_2": (ab["t_cut"], ab["arc"].Y[-1]),
                      "disposal": (dsp["t_ign"], dsp["arc"].Y[0])}.items():
        ev[k]["u_deg"] = float(C.elements(s[:3], s[3:6])["u"])
    prop_used = s2["prop_used"] + ab["prop"] + dsp["prop"]
    el = lambda e: {k: float(e[k]) for k in ("a", "e", "i", "raan", "argp", "u", "hp", "ha", "period")}  # noqa: E731
    fl = {k: {kk: f[kk] for kk in ("t_sep", "t", "lat", "lon", "range_km")} for k, f in sol["falls"].items()}
    for k, f in fl.items():
        f["zone"] = [z["id"] for z in sol["zones"] if C.point_in_ring(f["lat"], f["lon"], z["ring"])]
        f["beta"] = {"boosters": v.beta_b, "core": v.beta_c, "fairing": v.beta_f}[k]
    s = Satrec.twoline2rv(sol["tle"][-2], sol["tle"][-1])
    return {
        "t0": fr.t0.isoformat(), "pad": {"lat": C.PAD[0], "lon": C.PAD[1], "name": "Hainan commercial LC1 (GCAT site record; its LC1 row is 3.1 km south, at sea)"},
        "vehicle": C.vehicle_dict(v),
        "plane": {"inc": 50.0, "raan_deg": sol["raan"], "K_site_deg": sol["k_site"],
                  "gmst_t0_deg": float(np.degrees(fr.g0))},
        "steering": {"az_deg": lw["az"], "kick_deg": lw["kick"], "t_q": lw["t_q"], "b_boost": lw["b_boost"],
                     "b_core": lw["b_core"], "tan0_core": lw["tan0_core"],
                     "zone_fit_resid_km": lw["fit_resid_km"], "zone_aim_km": lw["aim_km"],
                     "s2_burn1": s2["x"].tolist(), "s2_tan0": s2["tan0"],
                     "s2_burn1_pitch_deg": [float(np.degrees(np.arctan(s2["x"][0]))),
                                            float(np.degrees(np.arctan(s2["x"][0] + s2["x"][1] * s2["x"][4])))],
                     "s2_burn2": {"tan_pitch": ab["tan_p"], "tan_yaw": ab["tan_y"]}},
        "max_q": {"q_kPa": float(q[iq] / 1e3), "t": float(lw["arc"].T[iq]), "mach": float(mach[iq])},
        "sightings": {"observer": {"lat": SG.OBS[0], "lon": SG.OBS[1], "h_km": SG.OBS_H}, "sigma_deg": SG.SIGMA_DEG,
                      "fit": lw["sightings"]},
        "events": ev,
        "falls": fl,
        "transfer_orbit": el(e_tr),
        "coast_s": ab["t_ign"] - s2["t_cut"],
        "apogee_burn": {"t_ign": ab["t_ign"], "t_cut": ab["t_cut"], "t_apo": ab["t_apo"], "t_burn": ab["t_burn"],
                        "dv_kms": float(np.linalg.norm(ab["arc"].Y[-1, 3:6]) - np.linalg.norm(ab["arc"].Y[0, 3:6]))},
        "parking_orbit_osculating": el(e_park),
        "parking_orbit_mean": {"alt_km": sol["mean_alt"], "inc": float(np.degrees(s.inclo)),
                               "raan": float(np.degrees(s.nodeo)), "ecc": float(s.ecco),
                               "argp": float(np.degrees(s.argpo)), "tle_fit_rms_km": sol["tle_rms"]},
        "stage_orbit_osculating": el(e_rb),
        "s2_budget": {"load": v.s2_prop, "burn1": s2["prop_used"], "burn2": ab["prop"], "disposal": dsp["prop"],
                      "used": prop_used, "left": v.s2_prop - prop_used, "burn1_s": s2["t_cut"] - v.t_s2,
                      "m_cut1": s2["m_cut"]},
        "tle": sol["tle"], "tle_rb": sol["tle_rb"],
        "zones": [{k: z[k] for k in ("id", "fir", "b", "c", "centroid", "range_km", "bearing")} for z in sol["zones"]],
        "profile": sol["profile"],
        "phase_shift_s": sol["shift"],
        "phase_lead_at_t0_deg": sol["lead"],
        "y12_reference": Y12,
        "assumed": {"T_deploy_after_burn2_s": T_DEPLOY,
                    "disposal": ("at the Y12 stage apogee argument of latitude" if sol["profile"] == "direct"
                                 else "300 s after deployment (propellant dump)" if sol["profile"] == "attached"
                                 else f"{T_DISPOSAL:.0f} s after deployment"),
                    "stage_perigee_km": RP_STAGE - C.RE,
                    "s2_sequence_s": ([C.T_CUT_DIRECT, sol["t_ign2"]] if sol["profile"] == "direct"
                                      else f"apogee on the inserted orbit at u = {sol['u_apo']:.1f} deg, vehicle "
                                           f"{sol['shift']:.1f} s ahead of the element set"
                                      if sol["profile"] == "attached"
                                      else f"cutoff at {sol['r_cut_alt']:.0f} km perigee"),
                    "target_mean_alt_km": Y12["mean_alt_km"]},
    }


def write_all(sol, tag="CZ-8A_Y13"):
    fr = sol["frame"]
    rows = []
    for ph, arc in segments(sol):
        step = 1.0 if arc.T[-1] <= 900 else 10.0
        for t, s in sample(arc, arc.T[0], arc.T[-1], step):
            if rows and float(rows[-1]["t"]) >= t - 1e-6:
                continue
            rows.append(row(fr, ph, t, s))
    write_csv(OUT / f"{tag}_traj.csv", rows)
    orows = [row(fr, "orbit", t, s) for t, s in sample(sol["sat"], sol["sat"].T[0], sol["sat"].T[-1], 10.0)]
    write_csv(OUT / f"{tag}_orbit.csv", orows)
    srows = []
    for arc, ph in ((sol["stage"], "stage coast"), (sol["dsp"]["arc"], "disposal"), (sol["rb"], "stage orbit")):
        srows += [row(fr, ph, t, s) for t, s in sample(arc, arc.T[0], arc.T[-1], 10.0 if ph != "disposal" else 1.0)]
    write_csv(OUT / f"{tag}_stage.csv", srows)
    frows = []
    for k, f in sol["falls"].items():
        frows += [row(fr, k, t, s) for t, s in sample(f["arc"], f["arc"].T[0], f["arc"].T[-1], 2.0)]
    write_csv(OUT / f"{tag}_falls.csv", frows)
    (OUT / f"{tag}_prelaunch.tle").write_text("\n".join(sol["tle"] + sol["tle_rb"]) + "\n")
    kml(sol, OUT / f"{tag}_traj.kml")
    js = summary(sol)
    (OUT / f"{tag}_solution.json").write_text(json.dumps(js, indent=1))
    return js


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t0", default=T0_NOMINAL.isoformat())
    ap.add_argument("--payload", type=float, default=C.VEH.payload)
    ap.add_argument("--profile", choices=("attached", "direct", "longcoast"), default="attached")
    ap.add_argument("--shift", type=float, default=None,
                    help="attached: seconds ahead of the inserted element set (default: solved from the S2 budget)")
    ap.add_argument("--u-apo", type=float, default=U_APO_ATTACHED, help="attached: apogee argument of latitude, deg")
    ap.add_argument("--full-thrust", action="store_true", help="the published 371 t vehicle, no core throttle")
    ap.add_argument("--sens", action="store_true")
    a = ap.parse_args()
    plan = json.loads((C.OUT / "j2_plan.json").read_text())
    t0 = datetime.fromisoformat(a.t0)
    veh = replace(C.VEH if a.full_thrust else VEH_NOM, payload=a.payload)
    print(f"T0 {t0} UTC, payload {veh.payload:.0f} kg, K {plan['K_y13_deg']:.3f}, profile {a.profile}")
    sol = solve(t0, veh, plan["K_y13_deg"], a.profile, u_apo=a.u_apo, shift=a.shift)
    js = write_all(sol)
    ev = js["events"]
    st = js["steering"]
    print(f"lower stages: lift-off {js['vehicle']['m0'] / 1e3:.1f} t, core {veh.core_thr:.3f} from T+{veh.t_thr:.0f} s | "
          f"az {st['az_deg']:.3f} kick {st['kick_deg']:.4f} t_q {st['t_q']:.1f} b_boost {st['b_boost']:+.3f} "
          f"b_core {st['b_core']:+.3f} | zone resid {np.round(st['zone_fit_resid_km'], 2)} km")
    for d in js["sightings"]["fit"]:
        print(f"  {d['label']:20s} seen az {d['az_obs']:5.1f} alt {d['alt_obs']:4.1f} | model T+{d['t']:6.1f} s "
              f"az {d['az']:6.2f} alt {d['alt']:5.2f} ({d['sep_deg']:.2f} deg off), {d['alt_km']:.1f} km")
    for k, f in js["falls"].items():
        print(f"  {k:9s} T+{f['t']:6.0f} s  {f['lat']:.3f} {f['lon']:.3f}  {f['range_km']:.0f} km  zone {f['zone']}")
    for k, e in ev.items():
        print(f"  {k:14s} T+{e['t']:7.1f}  {e['alt_km']:7.1f} km  v {e['v_inertial']:.3f}  fpa {e['fpa_deg']:6.2f}  "
              f"{e['lat']:7.3f} {e['lon']:8.3f}")
    tr, pk, mn = js["transfer_orbit"], js["parking_orbit_osculating"], js["parking_orbit_mean"]
    if js["phase_shift_s"] is not None:
        print(f"phase shift {js['phase_shift_s']:+.1f} s ahead of the inserted element set; S2 burn-1 pitch "
              f"{js['steering']['s2_burn1_pitch_deg'][0]:.0f} -> {js['steering']['s2_burn1_pitch_deg'][1]:.0f} deg")
    print(f"transfer {tr['hp']:.1f} x {tr['ha']:.1f} km i {tr['i']:.4f}; coast {js['coast_s']:.0f} s "
          f"({js['coast_s'] / 60:.2f} min); "
          f"burn2 {js['apogee_burn']['t_burn']:.1f} s dv {js['apogee_burn']['dv_kms'] * 1e3:.0f} m/s")
    print(f"parking osc {pk['hp']:.1f} x {pk['ha']:.1f}; mean alt {mn['alt_km']:.2f} i {mn['inc']:.4f} "
          f"raan {mn['raan']:.4f} e {mn['ecc']:.5f}")
    print("S2 budget", {k: round(v, 1) for k, v in js["s2_budget"].items()})
    print("max q", js["max_q"])
    print(f"phase lead at T0 over the pad {js['phase_lead_at_t0_deg']:+.2f} deg (Y12 catalog -16.83, Y10 +4.55)")
    print("\n".join(js["tle"] + js["tle_rb"]))
    if a.sens:
        sens = []
        K = plan["K_y13_deg"]
        for lab, vv, k, kw in (("payload +3 t", replace(veh, payload=veh.payload + 3000), K, {}),
                               ("relight 1950 s", veh, K, dict(t_ign2=1950.0)),
                               ("S2 176 kN, 23.5 t", replace(veh, s2_f=176e3, s2_prop=23500.0), K, {}),
                               ("beta core 200", replace(veh, beta_c=200.0), K, {}),
                               ("K Y10/Y11 (308.58)", veh, 308.584, {}),
                               ("long-coast profile (Y10 type)", veh, K, dict(profile="longcoast"))):
            kw.setdefault("profile", a.profile)
            s = solve(t0, vv, k, verbose=False, **kw)
            j = summary(s)
            sens.append(dict(case=lab, kick=j["steering"]["kick_deg"], az=j["steering"]["az_deg"],
                             b_core=j["steering"]["b_core"], t_cut1=j["events"]["s2_cutoff_1"]["t"],
                             cut_lat=j["events"]["s2_cutoff_1"]["lat"], cut_lon=j["events"]["s2_cutoff_1"]["lon"],
                             coast_s=j["coast_s"], t_burn2=j["apogee_burn"]["t_ign"],
                             burn2_lat=j["events"]["s2_ignition_2"]["lat"], burn2_lon=j["events"]["s2_ignition_2"]["lon"],
                             s2_left=j["s2_budget"]["left"], fairing_km=j["falls"]["fairing"]["range_km"],
                             fairing_zone=j["falls"]["fairing"]["zone"], max_q=j["max_q"]["q_kPa"],
                             lead=j["phase_lead_at_t0_deg"]))
            print(sens[-1])
        (OUT / "sensitivity.json").write_text(json.dumps(sens, indent=1))


if __name__ == "__main__":
    main()
