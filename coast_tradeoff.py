"""Attached long coast from the revision 5 lower stages: phase shift and S2 propellant left against the coast length,
the T0 second and the S2 load, plus the full-thrust (371 t) vehicle for comparison.

Every case flies the lower stages fitted to both NOTAM zones and the four ground sightings (cz8a.fit_lower,
sight=True), a lofted S2 burn 1, a coast to an apogee on the inserted orbit and burn 2 onto its velocity
(build_traj.attach).  The shift (s ahead of outside-ref/outside.tle) is 0 whenever the S2 load can pay for it,
otherwise the smallest one that leaves build_traj.S2_RESERVE after burn 2.

    ~/.venvs/astro313/bin/python coast_tradeoff.py   -> output/coast_tradeoff.json
"""
import json
from dataclasses import replace
from datetime import datetime, timedelta

import numpy as np

import build_traj as B
import cz8a as C

T0_MIN = datetime(2026, 11, 20, 11, 59, 0)
_LOWER = {}


def lower(fr, veh, key):
    """The lower-stage fit is Earth-fixed, so it does not depend on T0: fit once per vehicle, reuse the solution."""
    x0 = _LOWER.get(key, B.X0_LOWER)
    lw = C.fit_lower(fr, C.load_zones(), veh, x0=x0, aim=B.AIM_KM, sight=True)
    _LOWER[key] = lw["fit_x"]
    return lw


def case(label, sec=59.0, u_apo=B.U_APO_ATTACHED, load=None, full=False):
    veh = replace(C.VEH if full else B.VEH_NOM, **({} if load is None else dict(s2_prop=load)))
    fr = C.Frame(T0_MIN + timedelta(seconds=sec))
    zones = C.load_zones()
    lw = lower(fr, veh, "full" if full else "nom")
    s2, ab, _, shift, _ = B.attach(fr, lw, u_apo, veh)
    fl = C.falls(fr, lw, veh, s2)
    cut = C.state_summary(fr, s2["t_cut"], s2["arc"].Y[-1])
    e = C.elements(s2["arc"].r_end, s2["arc"].v_end)
    dsp_kg = 41.0                       # disposal dump, as in build_traj
    x = s2["x"]
    row = dict(case=label, t0_s=sec, u_apo_deg=u_apo, s2_load_kg=veh.s2_prop, full_thrust=full, liftoff_t=veh.m0() / 1e3,
               shift_s=shift, s2_left_kg=veh.s2_prop - s2["prop_used"] - ab["prop"] - dsp_kg,
               coast_min=(ab["t_ign"] - s2["t_cut"]) / 60.0, t_cut1=s2["t_cut"], cut_alt_km=cut["alt_km"],
               transfer_hp_km=float(e["hp"]), transfer_ha_km=float(e["ha"]), t_ign2=ab["t_ign"], burn2_s=ab["t_burn"],
               pitch_deg=[float(np.degrees(np.arctan(x[0]))), float(np.degrees(np.arctan(x[0] + x[1] * x[4])))],
               booster_sep_alt_km=C.state_summary(fr, veh.t_bsep, lw["y_bsep"])["alt_km"],
               fairing_alt_km=C.state_summary(fr, veh.t_fair, s2["arc"](veh.t_fair))["alt_km"],
               sightings_max_deg=max(d["sep_deg"] for d in lw["sightings"]),
               zone_resid_max_km=max(abs(v) for v in lw["fit_resid_km"]),
               falls={k: dict(t=f["t"], range_km=f["range_km"],
                              zone=[z["id"] for z in zones if C.point_in_ring(f["lat"], f["lon"], z["ring"])])
                      for k, f in fl.items()})
    print(f"{label:30s} shift {shift:+6.1f} s | S2 left {row['s2_left_kg']:5.0f} kg | coast {row['coast_min']:5.2f} min | "
          f"cutoff T+{s2['t_cut']:.0f} s {cut['alt_km']:.0f} km, transfer {e['hp']:.0f} x {e['ha']:.0f} | pitch "
          f"{row['pitch_deg'][0]:.0f} -> {row['pitch_deg'][1]:.0f} deg | sightings <= {row['sightings_max_deg']:.2f} deg | "
          f"lift-off {row['liftoff_t']:.1f} t", flush=True)
    return row


def main():
    rows = [case("nominal")]
    rows += [case(f"coast, apogee u {u:.0f}", u_apo=u) for u in (330.0, 337.0, 344.0, 351.0)]
    rows += [case(f"T0 11:59:{sec:02.0f}", sec=sec) for sec in (0.0, 30.0, 45.0)]
    rows += [case(f"S2 load {load / 1e3:.1f} t", load=load) for load in (27200.0, 28700.0)]
    rows += [case("full thrust (371 t)", full=True)]
    (C.OUT / "coast_tradeoff.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
