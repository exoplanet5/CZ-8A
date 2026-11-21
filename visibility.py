"""Sunlight along the CZ-8A Y13 flight and the 10 deg ground visibility zone of the powered ascent.

The launch is at 19:59 Beijing time, so the whole powered ascent runs in Earth's shadow (the Sun is ~40 deg below the
horizon under S2 cutoff, a shadow height near 1800 km).  No sunlit plume, so no twilight jelly: the zone here is the
self-luminous engine flame, i.e. the union over the powered ascent (lift-off .. S2 cutoff 1) of ground points that see
the vehicle at >= 10 deg elevation while their own Sun-centre altitude is <= 0 deg at that instant, clipped to the
display night side at the reference epoch (S2 cutoff 1).  Sunlight: geometric Sun centre against the WGS84 ellipsoid
(JPL DE440s from ~/JPLDE), the convention of ../launch_traj_6A_reverse/ascent_visibility.py.

    ~/.venvs/astro313/bin/python visibility.py   -> output/visibility.json, output/illumination.csv
"""
import csv
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from skyfield.api import Loader
from skyfield.framelib import itrs

import cz8a as C

ELEVATION_DEG = 10.0
EXTENT = [103.0, 148.0, -22.0, 27.0]         # lon west, east, lat south, north (= make_map.EXTENT)
F = 1 / 298.257223563
AXES = np.array([C.RE, C.RE, C.RE * (1 - F)])
OUT = C.OUT
load = Loader(str(Path("~/JPLDE").expanduser()))
ts = load.timescale(builtin=True)
eph = load("de440s.bsp")


def sun_ecef(epoch, elapsed):
    t = ts.from_datetime(epoch.replace(tzinfo=timezone.utc))
    when = ts.tt_jd(t.tt + np.asarray(elapsed, float) / 86400.0)
    return (eph["sun"] - eph["earth"]).at(when).frame_xyz(itrs).km.T


def clearance(pos, sun):
    """> 0 when the segment vehicle -> Sun centre misses the ellipsoid."""
    r = np.asarray(pos) / AXES
    d = (np.asarray(sun) - pos) / AXES
    f = np.clip(-np.sum(r * d, axis=-1) / np.sum(d * d, axis=-1), 0.0, 1.0)
    c = r + np.asarray(f)[..., None] * d
    return np.sum(c * c, axis=-1) - 1.0


def sun_alt(surface, zenith, sun):
    ray = sun - surface
    return np.degrees(np.arcsin(np.clip(np.sum(ray * zenith, axis=-1) / np.linalg.norm(ray, axis=-1), -1, 1)))


def read(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: (np.array([r[k] for r in rows]) if k in ("utc", "phase") else np.array([float(r[k]) for r in rows]))
            for k in rows[0]}


def grid(extent, step):
    w, e, s, n = extent
    lon, lat = np.meshgrid(np.linspace(w, e, int(round((e - w) / step)) + 1), np.linspace(s, n, int(round((n - s) / step)) + 1))
    surf = np.moveaxis(C.geodetic_to_ecef(lat, lon, 0.0), 0, -1)
    la, lo = np.radians(lat), np.radians(lon)
    zen = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)
    return lon, lat, surf, zen


def compute(extent=EXTENT, step=0.08, dt=2.0):
    sol = json.loads((OUT / "CZ-8A_Y13_solution.json").read_text())
    t0 = datetime.fromisoformat(sol["t0"])
    tr = read(OUT / "CZ-8A_Y13_traj.csv")
    t_cut = sol["events"]["s2_cutoff_1"]["t"]
    t_end = sol["apogee_burn"]["t_cut"]
    xyz = np.c_[tr["x"], tr["y"], tr["z"]]
    pos = lambda t: np.array([np.interp(t, tr["t"], xyz[:, k]) for k in range(3)]).T  # noqa: E731
    # sunlight along the whole flight to the end of the apogee burn
    ts_all = np.r_[np.arange(0.0, t_end, 5.0), t_end]
    cl = clearance(pos(ts_all), sun_ecef(t0, ts_all))
    k = np.flatnonzero((cl[:-1] <= 0) & (cl[1:] > 0))
    pickup = None
    if len(k):
        tp = brentq(lambda t: float(clearance(pos(t), sun_ecef(t0, t))), ts_all[k[0]], ts_all[k[0] + 1], xtol=1e-4)
        lat, lon, h = C.ecef_to_geodetic(pos(tp))
        pickup = dict(t_s=float(tp), utc=(t0 + timedelta(seconds=tp)).isoformat(timespec="seconds"),
                      lat=float(lat), lon=float(lon), alt_km=float(h),
                      phase=str(tr["phase"][np.searchsorted(tr["t"], tp)]))
    powered_sunlit = bool(np.any(cl[ts_all <= t_cut] > 0))
    # Sun depression under the vehicle at key events
    dep = {}
    for key in ("booster_sep", "s2_ignition", "s2_cutoff_1"):
        e = sol["events"][key]
        g = C.geodetic_to_ecef(e["lat"], e["lon"], 0.0)
        la, lo = np.radians(e["lat"]), np.radians(e["lon"])
        z = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
        dep[key] = float(sun_alt(g, z, sun_ecef(t0, e["t"])))
    # 10 deg zone of the self-luminous powered ascent
    lon, lat, surf, zen = grid(extent, step)
    st = np.r_[np.arange(0.0, t_cut, dt), t_cut]
    sp, ss = pos(st), sun_ecef(t0, st)
    best = np.full(lon.shape, -90.0)
    for p, s in zip(sp, ss):
        ray = p - surf
        el = np.degrees(np.arcsin(np.clip(np.sum(ray * zen, axis=-1) / np.linalg.norm(ray, axis=-1), -1, 1)))
        dark = sun_alt(surf, zen, s) <= 0.0
        best = np.maximum(best, np.where(dark, el, -90.0))
    ref_sun = sun_alt(surf, zen, sun_ecef(t0, t_cut))
    field = np.minimum(best - ELEVATION_DEG, -ref_sun)
    visible = field >= 0
    cell = (np.radians(step) * C.RE) ** 2 * np.cos(np.radians(lat))
    summary = dict(
        definition=("Union over the powered ascent (T+0 .. S2 cutoff 1) of sea-level points seeing the vehicle at "
                    f">= {ELEVATION_DEG:.0f} deg elevation with their own Sun-centre altitude <= 0 deg at the same "
                    "instant; clipped to the night side at the reference epoch. Self-luminous flame, not sunlight."),
        illumination_model="Geometric Sun centre vs WGS84 ellipsoid, JPL DE440s; no refraction, terrain or plume size",
        elevation_deg=ELEVATION_DEG, reference_epoch_utc=(t0 + timedelta(seconds=t_cut)).isoformat(timespec="seconds"),
        extent=extent, grid_deg=step, time_step_s=dt, powered_interval_s=[0.0, t_cut],
        powered_ascent_sunlit=powered_sunlit, sunlight_pickup=pickup, sun_alt_under_vehicle_deg=dep,
        zone_area_km2=float(cell[visible].sum()),
        zone_extent=dict(lon=[float(lon[visible].min()), float(lon[visible].max())],
                         lat=[float(lat[visible].min()), float(lat[visible].max())]),
        checks=dict(powered_ascent_in_shadow=not powered_sunlit, nonempty=bool(visible.any()),
                    zone_inside_extent=bool(not (visible[0].any() or visible[-1].any() or visible[:, 0].any()
                                                 or visible[:, -1].any())),
                    no_daylight_cell=bool(np.all(ref_sun[visible] <= 0))))
    (OUT / "visibility.json").write_text(json.dumps(summary, indent=1))
    with open(OUT / "illumination.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_s", "utc", "lat", "lon", "alt_km", "sun_centre_visible", "ellipsoid_ray_clearance"])
        for t, c in zip(ts_all, cl):
            la, lo, h = C.ecef_to_geodetic(pos(t))
            w.writerow([f"{t:.1f}", (t0 + timedelta(seconds=float(t))).isoformat(timespec="seconds"),
                        f"{la:.4f}", f"{lo:.4f}", f"{h:.2f}", int(c > 0), f"{c:.6f}"])
    return dict(lon=lon, lat=lat, field=field, ref_sun=ref_sun, summary=summary)


if __name__ == "__main__":
    r = compute()
    print(json.dumps(r["summary"], indent=1))
