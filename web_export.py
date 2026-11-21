#!/usr/bin/env python3
"""Export the CZ-8A Y13 model to the 3D page's mission data (schema: ~/launch/launch-page-template/PROFILE_FORMAT.md).

Resamples the built trajectory onto uniform grids; it fits and integrates nothing.  The Earth-fixed track holds for every
T0 on the page; only the Sun moves.  Sun and Moon vectors come from the shared DE440s kernel, which is read, never copied.

Run after build_traj.py and make_map.py (visibility.json gives the sunlight pickup):
    ~/.venvs/astro313/bin/python web_export.py                 # -> docs/data/cz8a-y13.json
    ~/.venvs/astro313/bin/python web_export.py --out FILE      # anywhere else (new_page.py copies it into a new page)
Then register or check it with the launch-page skill's new_page.py (--add for an existing page).
"""
import argparse
import csv
import json
import math
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from skyfield.api import Loader, load_file
from skyfield.framelib import itrs

import cz8a as C

STEM = "CZ-8A_Y13"
MISSION_ID = "cz8a-y13"
REPO = "https://github.com/exoplanet5/CZ-8A"
README = f"{REPO}/blob/main/missions/{C.MISSION.name}/README.md"
T0_MIN, T0_MAX = "11:51:00", "12:17:00"       # B6120/26 window (B6122/26 closes 12:28)
REVS = 2.0                                     # parking-orbit revolutions after burn 2
SIDEREAL_DEG_PER_MIN = 360.98564736629 / 1440.0


def read(name):
    with open(C.OUT / name) as fh:
        rows = list(csv.DictReader(fh))
    t = np.array([float(q["t"]) for q in rows])
    X = np.array([[float(q[k]) for k in ("x", "y", "z")] for q in rows])
    return t, X, [q["phase"] for q in rows]


def uniform(t, X, t_a, t_b, dt_max, cubic=False):
    """Uniform samples t_a..t_b (both exact) from the source rows inside that span, ECEF km -> geodetic lon/lat/alt.
    Linear for the 1 s powered rows; a cubic spline through the 10 s coast and orbit rows (chord error 0.1 km linear)."""
    tol = 2e-3                                             # the CSVs round times to 1 ms; event times are exact
    keep = (t >= t_a - tol) & (t <= t_b + tol)
    ts, Xs = t[keep], X[keep]
    ts, first = np.unique(ts, return_index=True)          # phase boundaries repeat a row with the new phase
    Xs = Xs[first]
    if abs(ts[0] - t_a) < tol:                             # snap boundary rows, so adjacent segments share them exactly
        ts[0] = t_a
    if abs(ts[-1] - t_b) < tol:
        ts[-1] = t_b
    n = max(1, int(math.ceil((t_b - t_a) / dt_max - 1e-9)))
    grid = np.linspace(t_a, t_b, n + 1)
    P = CubicSpline(ts, Xs, axis=0)(grid) if cubic and len(ts) > 3 else np.column_stack([np.interp(grid, ts, Xs[:, k]) for k in range(3)])
    lon, lat, alt = [], [], []
    for p in P:
        la, lo, h = C.ecef_to_geodetic(p)
        lat.append(round(float(la), 6))
        lon.append(round(float(lo), 6))
        alt.append(round(max(float(h), 0.0), 4))
    return dict(met0=float(t_a), dt=float((t_b - t_a) / n), lon=lon, lat=lat, alt=alt)


def fall_track(t, X, ph, key, color, label):
    m = np.array([p == key for p in ph])
    tr = uniform(t[m], X[m], t[m][0], t[m][-1], 5.0)
    tr.update(color=color, label=label, dashed=True)
    return tr


def body_table(body, start, end, model):
    load = Loader(str(Path("~/JPLDE").expanduser()), verbose=False)
    ts = load.timescale(builtin=True)
    eph = load_file(str(Path("~/JPLDE/de440s.bsp").expanduser()))
    n = int((end - start).total_seconds() // 30) + 2
    ta = ts.from_datetimes([start + timedelta(seconds=30 * k) for k in range(n)])
    vec = eph["earth"].at(ta).observe(eph[body]).frame_xyz(itrs).km.T
    return dict(start_utc=start.strftime("%Y-%m-%dT%H:%M:%SZ"), dt=30, vectors=[[round(float(v), 1) for v in p] for p in vec], model=model)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=C.HERE / "docs" / "data" / f"{MISSION_ID}.json")
    args = ap.parse_args()

    sol = json.loads((C.OUT / f"{STEM}_solution.json").read_text())
    vis = json.loads((C.OUT / "visibility.json").read_text())
    ev, fl, sg = sol["events"], sol["falls"], sol["sightings"]
    t0 = datetime.fromisoformat(sol["t0"]).replace(tzinfo=timezone.utc)
    t_bsep, t_ccut, t_csep = ev["booster_sep"]["t"], ev["core_cutoff"]["t"], sol["vehicle"]["t_csep"]
    t_s2, t_fair, t_cut1 = ev["s2_ignition"]["t"], ev["fairing"]["t"], ev["s2_cutoff_1"]["t"]
    t_ign2, t_cut2, t_dep, t_disp = ev["s2_ignition_2"]["t"], ev["s2_cutoff_2"]["t"], ev["deployment"]["t"], ev["disposal"]["t"]
    park, park_m, tro = sol["parking_orbit_osculating"], sol["parking_orbit_mean"], sol["transfer_orbit"]
    period = park["period"]
    t_end = t_cut2 + REVS * period

    # one time series: the flight to deployment (traj.csv), then the payload orbit (orbit.csv)
    ta, Xa, _ = read(f"{STEM}_traj.csv")
    to, Xo, _ = read(f"{STEM}_orbit.csv")
    assert to[-1] >= t_end, "orbit.csv is shorter than the exported revolutions"
    t, X = np.r_[ta, to], np.r_[Xa, Xo]
    o = np.argsort(t, kind="stable")
    t, X = t[o], X[o]
    tracks = dict(ascent=uniform(t, X, 0.0, t_cut1, 1.0), coast=uniform(t, X, t_cut1, t_ign2, 5.0, cubic=True),
                  burn2=uniform(t, X, t_ign2, t_cut2, 0.5), orbit=uniform(t, X, t_cut2, t_end, 5.0, cubic=True))
    tracks["ascent"].update(color="#d95926", label="Ascent")
    tracks["coast"].update(color="#3987e5", label="Long coast")
    tracks["burn2"].update(color="#199e70", label="Apogee burn")
    tracks["orbit"].update(color="#c77dff", label="Parking orbit")
    tf, Xf, pf = read(f"{STEM}_falls.csv")
    tracks.update(boosterfall=fall_track(tf, Xf, pf, "boosters", "#ff8a50", "Booster fall"),
                  corefall=fall_track(tf, Xf, pf, "core", "#ffb84f", "Core fall"),
                  fairfall=fall_track(tf, Xf, pf, "fairing", "#cfd6dd", "Fairing fall"))

    def zone(z, role, label, color, desc, impact):
        ring = [[round(lo, 5), round(la, 5)] for la, lo in z["ring"]]
        if ring[0] != ring[-1]:
            ring.append(ring[0])
        return dict(key=z["id"].split("/")[0], id=z["id"], role=role, label=label, color=color, t0=z["t0"], t1=z["t1"], ring=ring,
                    centroid=[round(z["centroid"][1], 4), round(z["centroid"][0], 4)], dist_km=round(z["range_km"], 1),
                    bearing=round(z["bearing"], 2), desc=desc, impact=dict(met=round(impact["t"], 1), lon=round(impact["lon"], 4), lat=round(impact["lat"], 4)))

    zs = {z["id"]: z for z in C.load_zones()}
    zones = [zone(zs["B6120/26"], "Boosters", "2 liquid boosters", "#ff8a50",
                  f"RPHI NOTAM, fall area of unburned debris. The model puts the two K2 boosters here (separation T+{t_bsep:.1f} s, fitted to the centroid).", fl["boosters"]),
             zone(zs["B6122/26"], "Core + fairing", "Core stage + fairing", "#ffb84f",
                  f"RPHI NOTAM. The core (separation T+{t_csep:.1f} s) is fitted to the centroid; the fairing halves (T+{t_fair:.1f} s, ballistic coefficient 30) land "
                  f"here at T+{fl['fairing']['t']:.0f} s without being fitted.", fl["core"])]

    def e(met, label, major=False, short=None, track=None):
        out = dict(met=round(float(met), 3), label=label, major=major, kind="ship")
        if short:
            out["short"] = short
        if track:
            out["track"] = track
        return out

    pick = vis["sunlight_pickup"]
    sep = {q["label"]: q for q in sg["fit"]}["stage-1 separation"]
    mq = sol["max_q"]
    events = [e(0, "Liftoff", True), e(mq["t"], f"Max-Q {mq['q_kPa']:.0f} kPa"),
              e(t_bsep, f"Booster separation, {ev['booster_sep']['alt_km']:.0f} km", True, "Boosters off"),
              e(t_ccut, "Core cutoff"), e(t_csep, f"Stage 1 separation, {sep['alt_km']:.0f} km (sighted from Lingshui)", True, "Stage 1 sep."),
              e(t_s2, "S2 ignition 1"), e(t_fair, f"Fairing jettison, {ev['fairing']['alt_km']:.0f} km"),
              e(t_cut1, f"S2 cutoff 1: {tro['hp']:.0f} × {tro['ha']:.0f} km transfer orbit", True, "S2 cutoff 1"),
              e(pick["t_s"], f"First sunlight on the vehicle, {pick['alt_km']:.0f} km"),
              e(t_ign2, f"S2 ignition 2, apogee burn ({sol['apogee_burn']['dv_kms'] * 1e3:.0f} m/s)", True, "Apogee burn"),
              e(t_cut2, f"S2 cutoff 2: {park['hp']:.0f} × {park['ha']:.0f} km parking orbit", True, "Parking orbit"),
              e(t_dep, "Satellite deployment (assumed +120 s)", True, "Deployment"),
              e(t_disp, "Stage disposal, propellant dump (assumed +300 s)"),
              e(t_cut2 + period, "1 revolution"), e(t_end, "2 revolutions"),
              e(fl["boosters"]["t"], "Booster impact, B6120/26", True, "Booster impact", "boosterfall"),
              e(fl["core"]["t"], "Core impact, B6122/26", True, "Core impact", "corefall"),
              e(fl["fairing"]["t"], "Fairing impact, B6122/26", False, None, "fairfall")]
    events.sort(key=lambda q: q["met"])

    veh, b = sol["vehicle"], sol["s2_budget"]
    st = sol["steering"]
    notes = [
        "Post-launch model, not telemetry and not an official ephemeris. CZ-8A Y13 (SatNet LEO group 30, nine satellites) lifted off from Hainan "
        "commercial pad LC1 at 19:59 Beijing time on 20 November 2026 (Xinhua). The model T0 is 11:59:59 UTC; the seconds are not published.",
        "Lower stages: published sequence (boosters off T+171.9 s, stage 1 separation T+198.5 s, S2 ignition T+200.8 s, fairing T+209.4 s), "
        "fitted to the centroids of RPHI NOTAMs B6120/26 (boosters) and B6122/26 (core and fairing) and to four ground sightings from "
        f"N 18°26′14″ E 110°04′28″ near Lingshui. The sightings fit to {max(q['sep_deg'] for q in sg['fit']):.2f}° at worst; the timed stage-1 "
        f"separation is seen at {sep['az_obs']}/{sep['alt_obs']}° and modelled at {sep['az']:.1f}/{sep['alt']:.1f}°.",
        f"Open point: the sightings and zones fit together only if the core runs at {veh['core_thr'] * 100:.1f} % thrust from T+{veh['t_thr']:.0f} s "
        f"to booster separation. Lift-off is then {veh['m0'] / 1e3:.1f} t against the published 371 t. The throttle stands in for whatever "
        "lowers the lower-stage energy; it is not claimed as fact.",
        f"S2 load {b['load'] / 1e3:.1f} t (assumed). Burn 1 lasts {b['burn1_s']:.0f} s to a {tro['hp']:.0f} × {tro['ha']:.0f} km transfer orbit, "
        f"followed by a {sol['coast_s'] / 60:.1f} min long coast and a {sol['apogee_burn']['t_burn']:.1f} s apogee burn on the inserted orbit; "
        f"{b['left'] / 1e3:.2f} t is left.",
        "The inserted orbit is CZ-8A Y12's (2026-243) first element set moved in plane to Y13. The model reaches it with zero phase shift; "
        f"its mean altitude is {park_m['alt_km']:.1f} km at i {park_m['inc']:.3f}°. Y12 flew the same profile. From the 1050 km parking orbit, "
        "J2 drift takes the group into the plane of 2025-295.",
        f"The page's perigee and apogee ({park['hp']:.1f} × {park['ha']:.1f} km) are osculating after burn 2; the TLEs are SGP4 mean elements fitted "
        "to the model orbit. Catalog numbers 99913/99914 and the designator 26999 are placeholders. Deployment +120 s and disposal +300 s after "
        "burn 2 are assumptions.",
        f"Sunlight: geometric WGS84 occultation of the Sun's centre with DE440s, 30 s cadence. The powered ascent is in Earth's shadow; "
        f"the vehicle first meets sunlight during the coast at T+{pick['t_s'] / 60:.1f} min over New Zealand, and the apogee burn is in daylight. "
        "Shifting T0 on this page keeps the Earth-fixed path and moves only the Sun.",
        "Pad 19.5970°N 110.9307°E (GCAT site record; GCAT's LC1 launch-point row lies 3.1 km south, at sea). "
        "Data: GCAT (J. McDowell, planet4589.org), Space-Track.org.",
    ]
    meta = dict(
        title="CZ-8A Y13 · SatNet 50°", vehicle="CZ-8A Y13", mission="SatNet LEO group 30", date=t0.strftime("%Y-%m-%d"),
        subtitle=f"SatNet LEO group 30 · i {park_m['inc']:.1f}° · {park_m['alt_km']:.0f} km · long coast",
        pad=[C.PAD[1], C.PAD[0]], pad_alt=C.PAD_H, pad_name="Hainan commercial LC1",
        t0_default=t0.strftime("%H:%M:%S"), t0_min=T0_MIN, t0_max=T0_MAX,
        t0_presets=[dict(label="Window 11:51", time=T0_MIN), dict(label="Lift-off 11:59:59", time=t0.strftime("%H:%M:%S")),
                    dict(label="Window end 12:17", time=T0_MAX)],
        met_default=200, insertion_met=t_cut2, period_s=period, orbit_revolutions=REVS,
        inc=round(park["i"], 4), perigee_km=round(park["hp"], 1), apogee_km=round(park["ha"], 1), az=round(st["az_deg"], 3),
        s4_cutoff_met=t_cut1, cutoff_label="S2 cutoff 1",
        sky_mark=dict(met=t_csep, label="Stage 1 sep."),
        model_badge="Post-launch model",
        about=["Fitted to the NOTAM drop zones, four ground sightings and Y12's inserted orbit; not telemetry.",
               "WGS84 geometry · JPL DE440s sunlight."],
        readme_url=README,
        tle_note="Model TLEs at the T0 above, retimed with it. Catalog numbers and designator are placeholders.",
        # Rocket cartoon. Core 3.35 m and K2 boosters 2.25 m (CZ-7 family), S2 12.375 m x 3.35 m, 5.2 m fairing; lengths approximate.
        stack=dict(total_m=50.5, burn_end=t_cut2,
                   segments=dict(core=dict(len=26.3, w=3.35, kind="liquid", band=True), s2=dict(len=12.4, w=3.35, kind="liquid"),
                                 fairing=dict(len=11.8, w=5.2, kind="fairing"), payload=dict(len=6.0, w=3.2, kind="payload")),
                   states=[dict(until=t_bsep, segs=["core", "s2", "fairing"], boosters=dict(count=2, len=26.3, w=2.25, burn="liquid"), zoom=1, burn="liquid"),
                           dict(until=t_csep, segs=["core", "s2", "fairing"], zoom=1, burn="liquid"),
                           dict(until=t_fair, segs=["s2", "fairing"], zoom=1.6, burn="liquid"),
                           dict(until=t_dep, segs=["s2", "payload"], zoom=2.2, burn="liquid"),
                           dict(segs=["payload"], zoom=3.2, burn=None)]),
        observer_home="Lingshui sighting",
        extra_cities=[["Lingshui sighting", 110.0744, 18.4372, 3],
                      ["Manila", 120.98, 14.60, 1], ["Ho Chi Minh City", 106.70, 10.78, 1], ["Da Nang", 108.22, 16.07, 2],
                      ["Puerto Princesa", 118.74, 9.74, 2], ["Kota Kinabalu", 116.07, 5.98, 2], ["Bandar Seri Begawan", 114.94, 4.89, 3],
                      ["Cebu", 123.89, 10.32, 2], ["Davao", 125.61, 7.07, 2], ["Zamboanga", 122.08, 6.91, 3], ["Manado", 124.84, 1.47, 2],
                      ["Ternate", 127.38, 0.79, 3], ["Ambon", 128.18, -3.70, 2], ["Jayapura", 140.72, -2.53, 2], ["Koror", 134.48, 7.34, 3],
                      ["Hagåtña", 144.75, 13.47, 2], ["Darwin", 130.84, -12.46, 1], ["Port Moresby", 147.18, -9.44, 1], ["Cairns", 145.77, -16.92, 2],
                      ["Brisbane", 153.03, -27.47, 1], ["Sydney", 151.21, -33.87, 1], ["Nouméa", 166.46, -22.27, 2], ["Auckland", 174.76, -36.85, 1],
                      ["Wellington", 174.78, -41.29, 1], ["Christchurch", 172.64, -43.53, 2], ["Hanga Roa", -109.43, -27.15, 2],
                      ["Santiago", -70.67, -33.45, 1], ["Lima", -77.04, -12.05, 1], ["Singapore", 103.82, 1.35, 1], ["Kuala Lumpur", 101.69, 3.14, 1],
                      ["Bangkok", 100.50, 13.76, 1], ["Jakarta", 106.85, -6.21, 1]],
        flight_segments=["ascent", "coast", "burn2", "orbit"],
        powered_intervals=[[0.0, t_ccut], [t_s2, t_cut1], [t_ign2, t_cut2]],
        phases=[dict(until=0, label="Pre-launch"), dict(until=t_bsep, label="Boosters + core"), dict(until=t_ccut, label="Core stage"),
                dict(until=t_s2, label="Staging"), dict(until=t_cut1, label="S2 burn 1"), dict(until=t_ign2, label="Long coast"),
                dict(until=t_cut2, label="Apogee burn"), dict(until=t_dep, label="Parking orbit, before deployment"),
                dict(until=t_end + 1, label="Parking orbit, 1050 km")],
        views=[["Global", 150, -12, 3.6], ["Hainan", 111.15, 18.9, 1.1], ["Ascent corridor", 115.5, 13.5, 1.45],
               ["Long coast", 160, -25, 2.4], ["Apogee burn", -90, -27, 1.6]],
        visibility=dict(extent=[103, 148, -22, 27], grid_step_deg=0.25, time_step_s=2, levels=[10]),
        notes=notes)
    lo = datetime.combine(t0.date(), datetime.strptime(T0_MIN, "%H:%M:%S").time(), timezone.utc)
    hi = datetime.combine(t0.date(), datetime.strptime(T0_MAX, "%H:%M:%S").time(), timezone.utc)
    t_a, t_b = lo - timedelta(minutes=10), hi + timedelta(seconds=t_end + 60)
    sun = body_table("sun", t_a, t_b, "DE440s geometric Sun centre; WGS84 occultation")
    moon = body_table("moon", t_a, t_b, "DE440s geometric Moon centre")
    tle = sol["tle"] + sol["tle_rb"]
    data = dict(schema_version=1, meta=meta, tracks=tracks, zones=zones, events=events,
                tle_ref=dict(t0=t0.strftime("%Y-%m-%dT%H:%M:%SZ"), raan_rate_deg_per_min=SIDEREAL_DEG_PER_MIN,
                             sets=[dict(name=tle[0], l1=tle[1], l2=tle[2]), dict(name=tle[3], l1=tle[4], l2=tle[5])]),
                sun=sun, moon=moon)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(data, separators=(",", ":"), ensure_ascii=False))
    print(f"wrote {args.out} ({args.out.stat().st_size / 1e3:.0f} kB): "
          + ", ".join(f"{k} {len(v['lon'])} pts dt {v['dt']:.2f}" for k, v in tracks.items()))


if __name__ == "__main__":
    main()
