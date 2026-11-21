"""CZ-8A Y13 (Guowang / WHDW, i = 50 deg) from Hainan commercial pad LC1: 3-DOF ascent, long coast, apogee burn.

Flight law
  lift-off .. booster separation   4 x YF-100, vertical rise, pitch kick, zero incidence on an Earth-fixed great circle
                                   to t_q, then linear-tangent pitch from the flight path (rate b_boost); optional core
                                   throttle (core_thr from t_thr to booster separation)
  booster sep .. core cutoff       2 x YF-100, the same linear-tangent law continued at rate b_core, no attitude step, no
                                   yaw (without t_q: zero incidence to booster separation, core pitch from the flight path)
  S2 burn 1                        2 x YF-75H, linear-tangent pitch and yaw about the target plane, to a transfer orbit
                                   whose apogee lies on the inserted orbit (solve_attached, nominal), or to a perigee
                                   cutoff (solve_s2) or a lofted sub-orbital arc (solve_direct)
  long coast                       40 min to the transfer apogee, J2-J4
  S2 burn 2 (apogee)               onto the inserted orbit's velocity (or circularisation in the target plane)
  deployment, disposal burn        retrograde burn lowering the stage perigee to the Y12 stage value (904 km)

Event times are the published CZ-8 / CZ-8A sequence (booster separation 171.9 s from CZ-8 Y1; stage separation 198.5 s,
S2 ignition 200.8 s and fairing 209.4 s from the CZ-8A sequence). Masses are consumption over that sequence at the
published YF-100 flow plus GCAT dry masses; assumed items are marked ASSUMED and recorded in the solution JSON.

Generic pieces are imported, never copied: frames, J2-J4 gravity, atmosphere, coast / fall, elements and the SGP4 fit
from ../launch_traj_NROL-97/nrol97.py; the guided-burn kernel, pad geometry and NOTAM parser from
../launch_traj_CERES-2/ceres2.py.

Units: km, km/s, s, kg, N. Inertial frame = TEME of date (GMST rotation of ECEF).
"""
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.append(str(Path("~/launch/launch_traj_NROL-97").expanduser()))   # appended: these folders hold their own
sys.path.append(str(Path("~/launch/launch_traj_CERES-2").expanduser()))    # build_traj.py / make_map.py
from nrol97 import (MU, RE, WE, Arc, Frame, alt_of, atmosphere, bearing, coast, ecef_to_geodetic, elements,  # noqa: E402
                    fall, gc_dist, geodetic_to_ecef, plane_normal, rotz, v_rel)
from ceres2 import (M_FREE, M_GUIDED, M_RISE, M_TURN, P_AE, P_AREA, SHAPE_LIQ, _guided_P, _motor_P, _y,  # noqa: E402
                    centroid, cross_track, march, pad_frame, parse_notams, point_in_ring)
# ceres2 puts NROL-97 back at the front of sys.path, where its make_map.py would shadow this project's: move it last
_NROL = str(Path("~/launch/launch_traj_NROL-97").expanduser())
sys.path[:] = [q for q in sys.path if q != _NROL] + [_NROL]

G0 = 9.80665
PAD = (19.5970, 110.9307)          # Hainan commercial pad 1: GCAT's site record (19.5968, 110.9308) and the surveyed pad
                                   # of ../launch_traj_CZ-10B.  GCAT's LC1 row (19.5685) lies 3.1 km south, 2.4 km at sea
PAD_H = 0.01                       # km
# Each flight keeps its inputs and results in missions/<flight>/; Y13 is the only one modelled so far
MISSION = HERE / "missions" / "Y13_20261120"
OUT = MISSION / "output"
NOTAM = MISSION / "NOTAM-CZ-8A-Y13"
# thrust windows run 1 us past the burn end, so the last RK4 stage never rounds to tau > 1 (thrust off), which made
# burn-time derivatives noisy and stalled the solvers
EPS_T = 1e-6


@dataclass(frozen=True)
class Vehicle:
    # YF-100 (GCAT / Wikipedia): 1200 kN SL, 1340 kN vac, Isp 300 / 335 s
    f_sl: float = 1200e3
    f_vac: float = 1340e3
    isp_vac: float = 335.0
    # K2-1 booster x2: GCAT dry 6.0 t; residual ASSUMED
    b_dry: float = 6000.0
    b_res: float = 1000.0
    # K3 core: GCAT dry 12.0 t; residual ASSUMED
    c_dry: float = 12000.0
    c_res: float = 1500.0
    # S2: 2 x YF-75H, 200 kN, Isp 442.6 s (Wikipedia/GCAT); dry incl. multi-function cabin and dispenser ASSUMED
    s2_f: float = 200e3
    s2_isp: float = 442.6
    s2_dry: float = 3600.0
    s2_prop: float = 30200.0      # ASSUMED load (3.35 m stage; the CZ-8A sequence's 508.6 s first burn alone uses 23.4 t);
                                  # 30.2 t since revision 5 (user, 2026-11-21): lift-off 371.0 t at full thrust, the
                                  # published mass.  26.0 t in revisions 1-2, 27.2 t in 3-4
    fairing: float = 2900.0       # 5.2 m fairing, ASSUMED
    payload: float = 6300.0       # 9 x 700 kg, GCAT estimate
    # published sequence, s after lift-off
    t_rise: float = 10.0          # vertical rise, ASSUMED
    # core throttle: the two core YF-100 run at core_thr from t_thr to booster separation (CZ-8: 77.5 % possible).
    # Propellant is consumption, so the throttle lowers the core load and the lift-off mass.  1.0 = full thrust.
    core_thr: float = 1.0
    t_thr: float = 10.0
    t_bsep: float = 171.9         # booster cutoff = separation (CZ-8 Y1 actual 171.8783 s)
    t_ccut: float = 196.5         # core cutoff, ASSUMED 2 s before separation
    t_csep: float = 198.5         # CZ-8A stage separation
    t_s2: float = 200.8           # CZ-8A S2 first ignition
    t_fair: float = 209.4         # CZ-8A fairing jettison
    # aerodynamics (Cd from the shared kernel's curve): fairing 5.2 m; boosters 2.25 m
    area_b: float = 21.24 + 2 * 3.976
    area_c: float = 21.24
    # ballistic coefficients, kg/m^2, ASSUMED
    beta_b: float = 300.0
    beta_c: float = 400.0
    beta_f: float = 30.0
    s2_blend: float = 10.0

    @property
    def mdot1(self):
        return self.f_vac / (self.isp_vac * G0)

    @property
    def ae1(self):
        return (self.f_vac - self.f_sl) / 101325.0

    @property
    def b_prop(self):
        return self.mdot1 * self.t_bsep

    @property
    def c_prop(self):
        """core propellant = consumption: full flow to core cutoff, less the throttled share before booster sep."""
        return 2 * self.mdot1 * (self.t_ccut - (1 - self.core_thr) * max(self.t_bsep - self.t_thr, 0.0))

    @property
    def mdot2(self):
        return self.s2_f / (self.s2_isp * G0)

    def m_upper(self):
        return self.s2_dry + self.s2_prop + self.fairing + self.payload

    def m0(self):
        return 2 * (self.b_dry + self.b_res + self.b_prop) + self.c_dry + self.c_res + self.c_prop + self.m_upper()


VEH = Vehicle()


# ---------------------------------------------------------------- zones
def load_zones(path=NOTAM):
    z = parse_notams(path)
    for d in z:
        d["centroid"] = centroid(d["ring"])
        d["range_km"] = float(gc_dist(PAD, d["centroid"]))
        d["bearing"] = float(bearing(PAD, d["centroid"]))
    return sorted(z, key=lambda d: d["range_km"])


# ---------------------------------------------------------------- lower stages
def _unit(v):
    return v / np.linalg.norm(v)


def _tan_fp(y):
    """tan of the air-relative flight-path angle of state y"""
    vra = v_rel(y[:3], y[3:6])
    u = y[:3] / np.linalg.norm(y[:3])
    return (vra @ u) / np.linalg.norm(vra - (vra @ u) * u)


def fly_lower(frame, az, kick, b_core, veh=VEH, dt=0.2, t_q=None, b_boost=0.0):
    """Lift-off -> core cutoff.  az: Earth-fixed great-circle azimuth at the pad; kick: pitch-over from the vertical
    (deg) at t_rise; zero incidence on that great circle to t_q.  From t_q a linear-tangent pitch starts at the flight
    path with rate b_boost (tan per 100 s) and continues through the core-only phase at rate b_core, with no attitude
    step at booster separation.  t_q = None: zero incidence to booster separation and a core pitch starting from the
    flight path there (the law of revisions 1-4).  The core engines run at veh.core_thr from veh.t_thr to booster
    separation."""
    up, east, north = pad_frame(PAD)
    hz = np.sin(np.radians(az)) * east + np.cos(np.radians(az)) * north
    ng = np.cross(up, hz)
    r0, v0 = frame.to_eci(geodetic_to_ecef(PAD[0], PAD[1], PAD_H), 0.0, np.zeros(3))
    y = _y(r0, v0, veh.m0())
    tq = veh.t_bsep if t_q is None else float(t_q)
    cuts = {0.0, veh.t_rise, tq, veh.t_bsep}
    if veh.core_thr < 1.0 and 0.0 < veh.t_thr < veh.t_bsep:
        cuts.add(veh.t_thr)
    cuts = sorted(c for c in cuts if c <= veh.t_bsep)
    arcs, nq, tan_q = [], None, None
    for a, b in zip(cuts[:-1], cuts[1:]):
        k = veh.core_thr if a >= veh.t_thr else 1.0
        mot = _Motor((2 + 2 * k) * veh.f_vac, (2 + 2 * k) * veh.mdot1, (2 + 2 * k) * veh.ae1, b - a + EPS_T)
        if b <= veh.t_rise:
            arc = march(a, y, b, dt, M_RISE, _motor_P(frame, mot, a, veh.area_b, pad=PAD, ng=ng), SHAPE_LIQ)
        elif b <= tq:
            y0 = y
            if a == veh.t_rise:                       # pitch kick, applied to the air-relative velocity
                vr = v_rel(y[:3], y[3:6])
                d = rotz(frame.theta(veh.t_rise)) @ (np.cos(np.radians(kick)) * up + np.sin(np.radians(kick)) * hz)
                y = y.copy()
                y[3:6] = y[3:6] - vr + np.linalg.norm(vr) * d
            arc = march(a, y, b, dt, M_TURN, _motor_P(frame, mot, a, veh.area_b, pad=PAD, ng=ng), SHAPE_LIQ)
            arc.Y[0] = y0
        else:
            if nq is None:                            # corridor plane frozen in inertial space at t_q
                nq, tan_q = _unit(rotz(frame.theta(tq)) @ ng), float(_tan_fp(y))
            P = _guided_P(frame, mot.f_mean, mot.mdot, a, b - a + EPS_T, nq,
                          np.array([tan_q + b_boost * (a - tq) / 100.0, b_boost, 0.0, 0.0]), area=veh.area_b)
            P[P_AE] = mot.ae
            arc = march(a, y, b, dt, M_GUIDED, P, SHAPE_LIQ)
        arcs.append(arc)
        y = arc.Y[-1].copy()
    y_bsep = y.copy()
    y[6] -= 2 * (veh.b_dry + veh.b_res)
    if nq is None:                                    # core alone: pitch from the air-relative flight path
        nh, tan0 = _unit(rotz(frame.theta(veh.t_bsep)) @ ng), float(_tan_fp(y))
    else:                                             # the booster-phase pitch continued
        nh, tan0 = nq, tan_q + b_boost * (veh.t_bsep - tq) / 100.0
    Pc = _guided_P(frame, 2 * veh.f_vac, 2 * veh.mdot1, veh.t_bsep, veh.t_ccut - veh.t_bsep + EPS_T, nh,
                   np.array([tan0, b_core, 0.0, 0.0]), area=veh.area_c)
    Pc[P_AE] = 2 * veh.ae1                   # back-pressure loss of the two core nozzles
    a3 = march(veh.t_bsep, y, veh.t_ccut, dt, M_GUIDED, Pc, SHAPE_LIQ)
    return dict(arc=Arc.join(arcs + [a3]), y_bsep=y_bsep, ng=ng, nh_core=nh, tan0_core=float(tan0), az=float(az),
                kick=float(kick), b_core=float(b_core), t_q=None if t_q is None else float(t_q), b_boost=float(b_boost))


@dataclass(frozen=True)
class _Motor:
    f_mean: float
    mdot: float
    ae: float
    t_burn: float = 1.0


def coast_to(y, t0, t1, dt=0.25):
    return march(t0, y, t1, dt, M_FREE, np.zeros(22))


def falls(frame, lw, veh=VEH, s2=None):
    """Impact points of the boosters, the core and (if the S2 arc is given) the fairing halves."""
    out = {}
    yb = lw["y_bsep"]
    f = fall(yb[:3], yb[3:6], veh.t_bsep, veh.beta_b, dt_hi=1.0, dt_lo=0.25)
    out["boosters"] = _impact(frame, f, veh.t_bsep)
    c = coast_to(lw["arc"].Y[-1], veh.t_ccut, veh.t_csep)
    f = fall(c.r_end, c.v_end, veh.t_csep, veh.beta_c, dt_hi=1.0, dt_lo=0.25)
    out["core"] = _impact(frame, f, veh.t_csep)
    if s2 is not None:
        st = s2["arc"](veh.t_fair)
        f = fall(st[:3], st[3:6], veh.t_fair, veh.beta_f, dt_hi=1.0, dt_lo=0.25)
        out["fairing"] = _impact(frame, f, veh.t_fair)
    return out


def _impact(frame, f, t_sep):
    lat, lon, _ = frame.lla(f.r_end, f.t_end)
    return dict(arc=f, t_sep=float(t_sep), t=float(f.t_end), lat=float(lat), lon=float(lon),
                range_km=float(gc_dist(PAD, (lat, lon))))


def fit_lower(frame, zones, veh=VEH, x0=(3.0, -0.30, 140.0), aim=(0.0, 0.0), sight=False):
    """Kick, core pitch rate and azimuth putting the boosters on the near zone centroid and the core on the far one.
    Residuals: along- and cross-track miss (km) of both, 4 equations for 3 unknowns.
    aim: along-track offset (km, downrange positive) of the booster and core aim points from the two centroids.
    sight: also fit the ground sightings (sightings.py, 5 residuals, 0.3 deg weight against 5 km for the zones) with
    two more unknowns, the guided-pitch start t_q and its rate b_boost; x0 = (kick, b_core, az, t_q, b_boost)."""
    import sightings as SG
    zb, zc = zones[0]["centroid"], zones[1]["centroid"]

    def build(p, dt):
        if sight:
            return fly_lower(frame, p[2], p[0], p[1], veh, dt=dt, t_q=p[3], b_boost=p[4])
        return fly_lower(frame, p[2], p[0], p[1], veh, dt=dt)

    def evaluate(lw):
        fl = falls(frame, lw, veh)
        zr = []
        for key, z, off in (("boosters", zb, aim[0]), ("core", zc, aim[1])):
            al_t, cr_t = cross_track(PAD, lw["az"], z)
            al, cr = cross_track(PAD, lw["az"], (fl[key]["lat"], fl[key]["lon"]))
            zr += [al - al_t - off, cr - cr_t]
        if not sight:
            return zr, [], []
        c = coast_to(lw["arc"].Y[-1], veh.t_ccut, veh.t_s2)
        sr, det = SG.residuals(frame, Arc.join([lw["arc"], c]), {"t_csep": veh.t_csep})
        return zr, sr, det

    def res(p):
        zr, sr, _ = evaluate(build(p, 0.25))
        return np.r_[np.array(zr) / 5.0, sr]

    if sight:
        lb = [0.02, -10.0, 130.0, veh.t_rise + 20.0, -5.0]
        ub = [1.5, 10.0, 150.0, veh.t_bsep - 5.0, 5.0]
        x0 = np.clip(np.array(x0, float), np.array(lb) + 1e-6, np.array(ub) - 1e-6)
        sol = least_squares(res, x0, diff_step=1e-4, x_scale=[0.05, 0.1, 0.5, 10.0, 0.1], bounds=(lb, ub),
                            max_nfev=300)
    else:
        sol = least_squares(res, np.array(x0, float), diff_step=1e-4, x_scale=[0.5, 0.05, 0.5], max_nfev=200)
    lw = build(sol.x, 0.2)
    zr, sr, det = evaluate(lw)
    lw["fit_resid_km"] = [float(v) for v in zr]
    lw["aim_km"] = [float(aim[0]), float(aim[1])]
    lw["sightings"] = det
    lw["fit_x"] = [float(v) for v in sol.x]
    return lw


# ---------------------------------------------------------------- S2
def s2_ignition_state(lw, veh=VEH):
    c = coast_to(lw["arc"].Y[-1], veh.t_ccut, veh.t_s2)
    y = c.Y[-1].copy()
    y[6] = veh.m_upper()
    return c, y


def fly_s2(frame, lw, nh, x, veh=VEH, dt=0.5):
    """S2 burn 1: x = [tan pitch0, tan pitch1, tan yaw0, tan yaw1, t_burn/100] about plane nh, blended from the flight
    path over s2_blend seconds.  Fairing drops at t_fair."""
    c, y = s2_ignition_state(lw, veh)
    u = y[:3] / np.linalg.norm(y[:3])
    vra = v_rel(y[:3], y[3:6])
    tan0 = (vra @ u) / np.linalg.norm(vra - (vra @ u) * u)
    t_b = x[4] * 100.0
    P = _guided_P(frame, veh.s2_f, veh.mdot2, veh.t_s2, t_b + EPS_T, nh, np.asarray(x[:4]), area=veh.area_c, tan0=tan0,
                  blend=veh.s2_blend)
    b1 = march(veh.t_s2, y, veh.t_fair, dt / 5, M_GUIDED, P, SHAPE_LIQ)
    y = b1.Y[-1].copy()
    y[6] -= veh.fairing
    P2 = P.copy()
    P2[P_AREA] = 0.0                          # no drag after the fairing goes
    b2 = march(veh.t_fair, y, veh.t_s2 + t_b, dt, M_GUIDED, P2, SHAPE_LIQ)
    arc = Arc.join([c, b1, b2])
    return dict(arc=arc, t_ign=veh.t_s2, t_cut=veh.t_s2 + t_b, m_cut=float(b2.Y[-1, 6]), x=np.asarray(x, float),
                prop_used=veh.mdot2 * t_b, tan0=float(tan0))


def solve_s2(frame, lw, nh, r_cut, r_apo, veh=VEH, x0=None):
    """Burn 1 ending at radius r_cut with zero flight-path angle (perigee) and apogee radius r_apo, in plane nh."""
    def res(p):
        a = fly_s2(frame, lw, nh, p, veh, dt=1.0)["arc"]
        r, v = a.r_end, a.v_end
        e = elements(r, v)
        rn = np.linalg.norm(r)
        return [(rn - r_cut) / 2.0, (r @ v) / rn * 100.0, (e["ra"] - r_apo) / 5.0, (r @ nh) / 2.0, (v @ nh) * 500.0]

    starts = [] if x0 is None else [np.asarray(x0, float)]
    starts += [np.array([0.3, -0.25, 0.0, 0.0, 4.5]), np.array([0.5, -0.4, 0.05, -0.02, 4.7]),
               np.array([0.2, -0.15, -0.05, 0.02, 4.3]), np.array([0.8, -0.6, 0.0, 0.0, 4.9])]
    best = None
    for s in starts:
        try:
            so = least_squares(res, s, diff_step=1e-5, max_nfev=200, ftol=1e-12, xtol=1e-12, gtol=1e-12)
        except (ValueError, FloatingPointError, ZeroDivisionError):
            continue
        if best is None or so.cost < best.cost:
            best = so
        if so.cost < 1e-8:
            break
    s2 = fly_s2(frame, lw, nh, best.x, veh)
    s2["cost"] = float(best.cost)
    return s2


def time_to_apogee(r, v):
    """Two-body time from the state to the next apogee, s."""
    e = elements(r, v)
    nu = np.radians(e["nu"])
    E = 2 * np.arctan2(np.sqrt(1 - e["e"]) * np.sin(nu / 2), np.sqrt(1 + e["e"]) * np.cos(nu / 2))
    M = (E - e["e"] * np.sin(E)) % (2 * np.pi)
    n = np.sqrt(MU / e["a"] ** 3)
    return float(((np.pi - M) % (2 * np.pi)) / n)


def solve_s2_timed(frame, lw, nh, t_cut, r_apo, t_apo, veh=VEH, x0=None):
    """Burn 1 of fixed length ending at t_cut (published sequence) on an orbit whose apogee radius is r_apo and is
    reached at t_apo, in plane nh.  4 steering tangents for 4 conditions; the perigee falls out."""
    tb = (t_cut - veh.t_s2) / 100.0

    def res(p):
        a = fly_s2(frame, lw, nh, np.r_[p, tb], veh, dt=1.0)["arc"]
        r, v = a.r_end, a.v_end
        e = elements(r, v)
        return [(e["ra"] - r_apo) / 5.0, (t_cut + time_to_apogee(r, v) - t_apo) / 5.0, (r @ nh) / 2.0,
                (v @ nh) * 500.0]

    starts = [] if x0 is None else [np.asarray(x0, float)]
    starts += [np.array([0.6, -0.15, 0.16, -0.08]), np.array([0.8, -0.2, 0.1, -0.05]),
               np.array([1.0, -0.3, 0.16, -0.08]), np.array([0.4, -0.05, 0.1, -0.05])]
    best = None
    for s in starts:
        try:
            so = least_squares(res, s, diff_step=1e-5, max_nfev=300, ftol=1e-12, xtol=1e-12, gtol=1e-12)
        except (ValueError, FloatingPointError, ZeroDivisionError):
            continue
        if best is None or so.cost < best.cost:
            best = so
        if so.cost < 1e-8:
            break
    s2 = fly_s2(frame, lw, nh, np.r_[best.x, tb], veh)
    s2["cost"] = float(best.cost)
    return s2


T_CUT_DIRECT = 709.4        # CZ-8A typical sequence: S2 first cutoff
T_IGN2_DIRECT = 1761.5      # CZ-8A typical sequence: S2 second ignition


def solve_direct(frame, lw, nh, r_apo, veh=VEH, t_cut=T_CUT_DIRECT, t_ign2=T_IGN2_DIRECT):
    """Direct-ascent profile (Y1-Y8 and Y12 type): burn 1 of the published length on a lofted, sub-orbital arc whose
    apex (radius r_apo) is reached so that the circularisation burn, centred on the apex, lights at t_ign2.
    Several pitch starts; of the converged solutions the lofted branch (highest perigee) is kept, because the flat
    branch skims at 100-160 km for 500 s.  -> (s2, ab)"""
    best = None
    for start in ([1.5, -0.3, 0.16, -0.08], [1.0, -0.2, 0.16, -0.08], [2.0, -0.4, 0.16, -0.08],
                  [0.6, -0.15, 0.16, -0.08]):
        x0, t_target = start, t_ign2 + 15.0
        for _ in range(6):
            # two-body apex timing differs from the J2-integrated apex by several seconds: correct the target by
            # the integrated ignition error until the burn lights at t_ign2
            s2 = solve_s2_timed(frame, lw, nh, t_cut, r_apo, t_target, veh, x0=x0)
            x0 = s2["x"][:4]
            ab = apogee_burn(frame, s2, nh, veh)
            if abs(ab["t_ign"] - t_ign2) < 0.3:
                break
            t_target += t_ign2 - ab["t_ign"]
        if s2["cost"] < 1e-6:
            hp = elements(s2["arc"].r_end, s2["arc"].v_end)["hp"]
            if best is None or hp > best[0]:
                best = (hp, s2, ab)
    if best is None:
        raise RuntimeError("no direct-ascent solution")
    return best[1], best[2]


class Shifted:
    """An SGP4 orbit read `shift` seconds late: the same plane, altitude and shape, the object `shift` s ahead."""

    def __init__(self, sat, shift):
        self.sat, self.shift = sat, float(shift)

    def sgp4(self, jd, fr):
        return self.sat.sgp4(jd, fr + self.shift / 86400.0)


def inserted_orbit(path=MISSION / "outside-ref" / "outside.tle", shift=0.0):
    """The inserted orbit the user supplied (Y12's first 2026-243E TLE moved to Y13), as an SGP4 object.
    shift (s) moves the object ahead along that orbit; its plane and altitude are kept."""
    from sgp4.api import Satrec
    L = Path(path).read_text().splitlines()
    sat = Satrec.twoline2rv(L[1], L[2])
    return Shifted(sat, shift) if shift else sat


def sgp4_teme(sat, frame, t):
    from nrol97 import jd_of
    _, r, v = sat.sgp4(jd_of(frame.utc(t)), 0.0)
    return np.array(r), np.array(v)


def solve_attached(frame, lw, ins, u_apo, veh=VEH, x0=None, iters=8):
    """Long coast attached to the inserted orbit `ins` (SGP4): burn 1 ends on a transfer orbit whose apogee
    (radius r_apo) lies where `ins` passes argument of latitude u_apo, at the moment `ins` is there, in its plane.
    Unknowns: 4 steering tangents + burn time.  The two-body apogee prediction from the cutoff state is corrected
    by the integrated (J2) apogee time and radius and by the along-track miss from `ins` at burn-2 ignition; burn 2
    ends on the inserted orbit's velocity.  `ins` may be a Shifted orbit (the vehicle ahead of the element set).
    -> (s2, ab, t_apo_target)"""
    # time at which the inserted orbit crosses u_apo (first crossing after T+40 min)
    ts = np.arange(2400.0, 2400.0 + 7000.0, 10.0)
    us = np.array([elements(*sgp4_teme(ins, frame, t))["u"] for t in ts])
    d = (us - u_apo + 180) % 360 - 180
    k = int(np.flatnonzero((d[:-1] < 0) & (d[1:] >= 0))[0])
    t_apo = ts[k] - d[k] / (d[k + 1] - d[k]) * (ts[k + 1] - ts[k])
    nh1 = elements(*sgp4_teme(ins, frame, 780.0))["hhat"]      # plane near cutoff (J2 regresses it ~0.13 deg/52 min)
    nh2 = elements(*sgp4_teme(ins, frame, t_apo))["hhat"]
    r_apo = float(np.linalg.norm(sgp4_teme(ins, frame, t_apo)[0]))
    bias_t = bias_u = bias_r = 0.0
    x0 = np.asarray(x0 if x0 is not None else [5.0, -1.1, 0.04, -0.02, 5.8], float)

    def res(p):
        a = fly_s2(frame, lw, nh1, p, veh, dt=1.0)["arc"]
        r, v = a.r_end, a.v_end
        e = elements(r, v)
        return [(e["ra"] + bias_r - r_apo) / 2.0, (a.t_end + time_to_apogee(r, v) + bias_t - t_apo) / 5.0,
                ((e["argp"] + 180 + bias_u - u_apo + 180) % 360 - 180) * 20.0, (r @ nh1) / 2.0, (v @ nh1) * 500.0]

    for _ in range(iters):
        best = None
        for s in (x0, [5.0, -1.1, 0.04, -0.02, 5.8], [0.6, -0.2, 0.16, -0.08, 5.6], [1.0, -0.3, 0.1, -0.05, 5.6]):
            try:
                so = least_squares(res, s, diff_step=1e-5, max_nfev=300, ftol=1e-12, xtol=1e-12, gtol=1e-12)
            except (ValueError, FloatingPointError, ZeroDivisionError):
                continue
            if best is None or so.cost < best.cost:
                best = so
            if so.cost < 1e-8:
                break
        x0 = best.x
        s2 = fly_s2(frame, lw, nh1, x0, veh)
        s2["cost"] = float(best.cost)
        ab = apogee_burn(frame, s2, nh2, veh, match=lambda t: sgp4_teme(ins, frame, t))
        sa = ab["arc"].Y[0, :6]                      # state at burn-2 ignition, a few seconds before apogee
        ri, vi = sgp4_teme(ins, frame, ab["t_ign"])
        e2b = elements(s2["arc"].r_end, s2["arc"].v_end)
        pred_t = s2["t_cut"] + time_to_apogee(s2["arc"].r_end, s2["arc"].v_end)
        bias_t += ab["t_apo"] - (pred_t + bias_t)
        bias_u += (elements(sa[:3], sa[3:6])["u"] - elements(ri, vi)["u"] + 180) % 360 - 180
        bias_r += np.linalg.norm(sa[:3]) - (e2b["ra"] + bias_r)
        miss = np.linalg.norm(sa[:3] - ri)
        if abs(ab["t_apo"] - t_apo) < 0.5 and miss < 1.0:
            break
    return s2, ab, float(t_apo)


def lead_at_t0(sat, t0):
    """Phase lead (deg) at T0 of an SGP4 orbit over the pad: argument of latitude of the orbit propagated back to T0
    minus that of the pad's inertial position in the same plane.  A catalog TLE gives the flight's actual value."""
    from nrol97 import jd_of
    fr = Frame(t0)
    _, r, v = sat.sgp4(jd_of(t0), 0.0)
    el = elements(np.array(r), np.array(v))
    rp = fr.to_eci(geodetic_to_ecef(*PAD), 0.0)
    h = el["hhat"]
    nod = np.cross([0, 0, 1.0], h)
    nod /= np.linalg.norm(nod)
    m = np.cross(h, nod)
    return float((el["u"] - np.degrees(np.arctan2(rp @ m, rp @ nod)) + 180) % 360 - 180)


def fly_burn(frame, y, t_ign, t_b, nh, tan_p, tan_y, veh=VEH, dt=0.2):
    P = _guided_P(frame, veh.s2_f, veh.mdot2, t_ign, t_b + EPS_T, nh, np.array([tan_p, 0.0, tan_y, 0.0]))
    return march(t_ign, y, t_ign + t_b, dt, M_GUIDED, P, SHAPE_LIQ)


def apogee_burn(frame, s2, nh, veh=VEH, match=None):
    """Coast from burn-1 cutoff to the transfer apogee; circularise there with a burn centred on apogee.
    match: optional callable t -> (r, v) of a target orbit; the burn then ends on that orbit's radial, horizontal and
    normal velocity instead of a circular one."""
    y = s2["arc"].Y[-1].copy()
    c = march(s2["t_cut"], y, s2["t_cut"] + 4000.0, 2.0, M_FREE, np.zeros(22))
    rv = np.einsum("ij,ij->i", c.Y[:, :3], c.Y[:, 3:6])
    k = int(np.nonzero((rv[:-1] > 0) & (rv[1:] <= 0))[0][0])
    # refine apogee time on the Hermite interpolant
    ts = np.linspace(c.T[k], c.T[k + 1], 201)
    t_apo = ts[np.argmin([abs(c(t)[:3] @ c(t)[3:6]) for t in ts])]
    m = y[6]

    def state(t):
        s = c(t)
        return np.concatenate([s, [m]])

    def res(p):
        tan_p, tan_y, tb = p
        t_ign = t_apo - tb / 2
        b = fly_burn(frame, state(t_ign), t_ign, tb, nh, tan_p, tan_y, veh, dt=0.5)
        r, v = b.r_end, b.v_end
        rn = np.linalg.norm(r)
        vh = np.linalg.norm(np.cross(r, v)) / rn
        if match is not None:
            rt, vt = match(b.t_end)
            rtn = np.linalg.norm(rt)
            return [((r @ v) / rn - (rt @ vt) / rtn) * 1000.0,
                    (vh - np.linalg.norm(np.cross(rt, vt)) / rtn) * 1000.0, ((v - vt) @ nh) * 1000.0]
        return [(r @ v) / rn * 1000.0, (vh - np.sqrt(MU / rn)) * 1000.0, (v @ nh) * 1000.0]

    so = least_squares(res, [0.0, 0.0, 17.0], diff_step=1e-4, x_scale=[1e-3, 1e-3, 1.0], ftol=1e-14, xtol=1e-14,
                       gtol=1e-14)
    tan_p, tan_y, tb = so.x
    t_ign = t_apo - tb / 2
    pre = Arc(*_cut(c, s2["t_cut"], t_ign, m))
    b = fly_burn(frame, state(t_ign), t_ign, tb, nh, tan_p, tan_y, veh)
    return dict(coast=pre, arc=b, t_apo=float(t_apo), t_ign=float(t_ign), t_cut=float(t_ign + tb), t_burn=float(tb),
                tan_p=float(tan_p), tan_y=float(tan_y), prop=float(veh.mdot2 * tb), m_cut=float(b.Y[-1, 6]),
                resid=so.fun.tolist())


def _cut(arc, t0, t1, m):
    """Resample arc on [t0, t1] at <= 2 s."""
    n = max(int(np.ceil((t1 - t0) / 2.0)), 1)
    T = np.linspace(t0, t1, n + 1)
    Y = np.array([np.concatenate([arc(t), [m]]) for t in T])
    return T, Y


def disposal(frame, y, t_ign, nh, rp_target, veh=VEH):
    """Retrograde burn (plane normal flipped, so the guided frame points backwards) lowering perigee to rp_target."""
    def rp(tb):
        b = fly_burn(frame, y, t_ign, tb, -nh, 0.0, 0.0, veh, dt=0.05)
        return elements(b.r_end, b.v_end)["rp"]

    lo, hi = 0.0, 20.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if rp(mid) > rp_target:
            lo = mid
        else:
            hi = mid
    tb = 0.5 * (lo + hi)
    b = fly_burn(frame, y, t_ign, tb, -nh, 0.0, 0.0, veh, dt=0.05)
    return dict(arc=b, t_ign=float(t_ign), t_burn=float(tb), prop=float(veh.mdot2 * tb))


# ---------------------------------------------------------------- diagnostics
def q_profile(arc):
    q, mach = [], []
    for k in range(len(arc.T)):
        r, v = arc.Y[k, :3], arc.Y[k, 3:6]
        rho, p, cs = atmosphere(alt_of(r))
        vr = np.linalg.norm(v_rel(r, v)) * 1e3
        q.append(0.5 * rho * vr * vr)
        mach.append(vr / cs)
    return np.array(q), np.array(mach)


def state_summary(frame, t, s):
    r, v = s[:3], s[3:6]
    lat, lon, h = frame.lla(r, t)
    vr = v_rel(r, v)
    rn = np.linalg.norm(r)
    return dict(t=float(t), utc=frame.utc(t).isoformat(timespec="milliseconds"), lat=float(lat), lon=float(lon),
                alt_km=float(h), v_inertial=float(np.linalg.norm(v)), v_rel=float(np.linalg.norm(vr)),
                fpa_deg=float(np.degrees(np.arcsin((r @ vr) / (rn * np.linalg.norm(vr))))),
                range_km=float(gc_dist(PAD, (lat, lon))))


def plane_raan(raan_deg, inc=50.0):
    return plane_normal(inc, raan_deg)


def vehicle_dict(veh=VEH):
    d = asdict(veh)
    d.update(b_prop=veh.b_prop, c_prop=veh.c_prop, m_upper=veh.m_upper(), m0=veh.m0(), mdot1=veh.mdot1,
             mdot2=veh.mdot2, liftoff_tw=4 * veh.f_sl / (veh.m0() * G0))
    return d
