"""Can a 45-50 min pre-insertion coast connect T0 (11:59:00) to the inserted orbit (outside-ref/outside.tle)?

For the zone-fitted ascent, S2 burn 1 is solved to a transfer orbit with apogee at the parking radius, cutoff at a
chosen altitude and flight-path angle; the coast is integrated to apogee and circularised there.  The phase lead at T0
(cz8a.lead_at_t0) of the resulting orbit is compared with the inserted orbit's -18.39 deg, and split into
    cutoff term = u_cut - n*t_cut - u_pad     (set by the ascent)
    coast term  = du - n*dt over the coast    (>= 0 for any coast ending at a 1050 km apogee from a low perigee)
The miss is converted to the T0 shift that would absorb it.

    ~/.venvs/astro313/bin/python longcoast_check.py   -> output/longcoast_check.json
"""
import json
import sys, numpy as np
from datetime import datetime
from scipy.optimize import least_squares
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cz8a as C
from nrol97 import fit_tle
from sgp4.api import Satrec
t0 = datetime(2026, 11, 20, 11, 59); fr = C.Frame(t0)
lw = C.fly_lower(fr, 139.4127, 0.51281, 0.04583)
nh = C.plane_raan((np.degrees(fr.g0) + 309.696) % 360)
R_APO = C.RE + 1052.1

def solve(h_cut, fpa, x0):
    def res(p):
        a = C.fly_s2(fr, lw, nh, p, dt=1.0)['arc']; r, v = a.r_end, a.v_end; rn = np.linalg.norm(r)
        e = C.elements(r, v)
        return [(rn - C.RE - h_cut) / 2, ((r @ v) / rn / np.linalg.norm(v) - np.sin(np.radians(fpa))) * 300,
                (e['ra'] - R_APO) / 5, (r @ nh) / 2, (v @ nh) * 500]
    best = None
    for s in [x0, [0.5, -0.12, 0.16, -0.08, 4.35], [1.0, -0.25, 0.16, -0.08, 4.6], [0.3, -0.08, 0.16, -0.08, 4.2]]:
        so = least_squares(res, s, diff_step=1e-5, max_nfev=200, ftol=1e-12, xtol=1e-12, gtol=1e-12)
        if best is None or so.cost < best.cost: best = so
        if so.cost < 1e-8: break
    return C.fly_s2(fr, lw, nh, best.x), best

rows = []
print(' h_cut  fpa | t_cut  u_cut  hp    | coast   burn2 T+   u_b2  | lead (inserted -18.39) | cutoff term + coast term | S2 left')
x0 = [0.5191, -0.1208, 0.1608, -0.0776, 4.3491]
for h_cut, fpa in [(150, 0), (200, 0), (300, 0), (450, 0), (600, 0), (200, -10), (300, 10), (400, 20)]:
    s2, so = solve(h_cut, fpa, x0)
    if so.cost > 1e-6:
        print(f"{h_cut:5d} {fpa:4d} | no solution (cost {so.cost:.1e})"); continue
    x0 = so.x
    e = C.elements(s2['arc'].r_end, s2['arc'].v_end)
    ab = C.apogee_burn(fr, s2, nh)
    sat = C.march(ab['t_cut'], ab['arc'].Y[-1], ab['t_cut'] + 13000, 5.0, C.M_FREE, np.zeros(22))
    tle, _, _ = fit_tle(fr, sat, ab['t_cut'], ab['t_cut'] + 13000, 99999, '26999A', 'X')
    sat_ = Satrec.twoline2rv(tle[-2], tle[-1]); lead = C.lead_at_t0(sat_, t0)
    n = np.degrees(sat_.no_kozai)
    cut_term = (e['u'] - n * s2['t_cut'] / 60 - 153.09 + 180) % 360 - 180
    left = C.VEH.s2_prop - s2['prop_used'] - ab['prop']
    rows.append(dict(h_cut_km=h_cut, fpa_deg=fpa, t_cut=float(s2['t_cut']), u_cut=float(e['u']), hp=float(e['hp']),
                     coast_min=float((ab['t_ign'] - s2['t_cut']) / 60), t_burn2=float(ab['t_ign']), lead=float(lead),
                     miss_deg=float(lead + 18.39), t0_shift_needed_s=float((lead + 18.39) / n * 60),
                     cutoff_term=float(cut_term), coast_term=float(lead - cut_term), s2_left_kg=float(left)))
    print(f"{h_cut:5d} {fpa:4d} | {s2['t_cut']:5.0f} {e['u']:6.1f} {e['hp']:6.0f} | {(ab['t_ign']-s2['t_cut'])/60:5.1f} min {ab['t_ign']:6.0f} {C.elements(ab['arc'].r_end, ab['arc'].v_end)['u']:6.1f} | {lead:+7.2f} (miss {lead+18.39:+6.2f} deg) | {cut_term:+6.2f} {lead-cut_term:+6.2f} | {left:5.0f} kg", flush=True)
json.dump(rows, open(C.OUT / 'longcoast_check.json', 'w'), indent=1)
