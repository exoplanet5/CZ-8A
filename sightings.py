"""Ground sightings of the CZ-8A Y13 ascent (supplied by the user, 2026-11-21) and their residuals against a trajectory.

Observer N 18 26 14, E 110 04 28 (Lingshui, south-east Hainan), 157 km SW of the pad.  Topocentric azimuth / altitude:
three positions on the climb (time not given) and the stage-1 (core / S2) separation, which the model times at
Vehicle.t_csep.  Angles are geometric (refraction ~0.07 deg at 13.5 deg altitude is ignored).

Residuals: for an untimed sighting, the signed angular distance (deg) from the sightline to the trajectory's sky path,
perpendicular to that path at its closest point; for the timed one, the azimuth (x cos alt) and altitude differences
at that time.
"""
import numpy as np

from nrol97 import ecef_to_geodetic, geodetic_to_ecef

OBS = (18 + 26 / 60 + 14 / 3600, 110 + 4 / 60 + 28 / 3600)       # deg
OBS_H = 0.02                                                     # km, ASSUMED
SIGHTINGS = [  # (label, az, alt, event)
    ("climb 1", 45.5, 13.5, None),
    ("climb 2", 72.4, 25.5, None),
    ("climb 3", 89.4, 28.6, None),
    ("stage-1 separation", 94.6, 28.9, "t_csep"),
]
SIGMA_DEG = 0.3          # weight of one sighting, ASSUMED


def _enu(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    up = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    east = np.array([-np.sin(lo), np.cos(lo), 0.0])
    return up, east, np.cross(up, east)


_O = np.array(geodetic_to_ecef(OBS[0], OBS[1], OBS_H))
_UP, _E, _N = _enu(*OBS)


def unit(az, alt):
    a, e = np.radians(az), np.radians(alt)
    return np.cos(e) * (np.sin(a) * _E + np.cos(a) * _N) + np.sin(e) * _UP


def azalt(r_ecef):
    """Topocentric azimuth and altitude (deg) of ECEF points (..., 3) from the observer."""
    d = np.asarray(r_ecef) - _O
    d = d / np.linalg.norm(d, axis=-1, keepdims=True)
    return (np.degrees(np.arctan2(d @ _E, d @ _N)) % 360.0, np.degrees(np.arcsin(np.clip(d @ _UP, -1, 1))))


def los(r_ecef):
    d = np.asarray(r_ecef) - _O
    return d / np.linalg.norm(d, axis=-1, keepdims=True)


def residuals(frame, arc, events, t_max=None):
    """Residual vector (deg) of `arc` (inertial, Arc with T, Y) against SIGHTINGS.  events: {'t_csep': t, ...}.
    -> (res list, detail list)"""
    T = arc.T if t_max is None else arc.T[arc.T <= t_max]
    X = np.array([frame.to_ecef(arc.Y[k, :3], T[k]) for k in range(len(T))])
    L = los(X)
    res, det = [], []
    for lab, az, alt, ev in SIGHTINGS:
        u = unit(az, alt)
        if ev is None:
            c = L @ u
            k = int(np.clip(np.argmax(c), 1, len(T) - 2))
            # refine on a parabola through the three samples around the best one
            t = T[k - 1:k + 2]
            p = np.polyfit(t - t[1], c[k - 1:k + 2], 2)
            tk = t[1] - p[1] / (2 * p[0]) if p[0] < 0 else t[1]
            tk = float(np.clip(tk, t[0], t[2]))
            rk = frame.to_ecef(arc(tk)[:3], tk)
            lk = los(rk)
            # sky-path direction at tk, and the signed offset of the sightline perpendicular to it
            dt = 1.0
            l2 = los(frame.to_ecef(arc(tk + dt)[:3], tk + dt))
            tang = l2 - lk
            tang -= (tang @ lk) * lk
            tang /= np.linalg.norm(tang)
            perp = np.cross(lk, tang)
            off = float(np.degrees(np.arcsin(np.clip(u @ perp, -1, 1))))
            res.append(off / SIGMA_DEG)
            a, e = azalt(rk)
            det.append(dict(label=lab, az_obs=az, alt_obs=alt, t=tk, az=float(a), alt=float(e), off_deg=off,
                            sep_deg=float(np.degrees(np.arccos(np.clip(lk @ u, -1, 1)))),
                            alt_km=float(ecef_to_geodetic(rk)[2]), range_km=float(np.linalg.norm(rk - _O))))
        else:
            tk = events[ev]
            rk = frame.to_ecef(arc(tk)[:3], tk)
            a, e = azalt(rk)
            daz = ((a - az + 180) % 360 - 180) * np.cos(np.radians(alt))
            res += [daz / SIGMA_DEG, (e - alt) / SIGMA_DEG]
            det.append(dict(label=lab, az_obs=az, alt_obs=alt, t=tk, az=float(a), alt=float(e), daz_cos=float(daz),
                            dalt=float(e - alt), sep_deg=float(np.degrees(np.arccos(np.clip(los(rk) @ u, -1, 1)))),
                            alt_km=float(ecef_to_geodetic(rk)[2]), range_km=float(np.linalg.norm(rk - _O))))
    return res, det
